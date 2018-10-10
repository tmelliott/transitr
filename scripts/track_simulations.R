source("scripts/common.R")

library(RSQLite)
library(dbplyr)

simnames <- list.files("simulations", pattern = "sim*", include.dirs = TRUE)
if (file.exists("simulations/arrivaldata.rda")) {
    load("simulations/arrivaldata.rda")
} else {
    tmax <- function(x, t) x[which.max(t)]
    pbapply::pboptions(type = "timer")
    arrivaldata <- pbapply::pblapply(list.files("simulations/archive", pattern = "trip_*", full.names = TRUE),
        function(f) {
            ent <- read(transit_realtime.FeedMessage, f)$entity
            do.call(bind_rows, 
                lapply(ent, function(e) {
                    stu <- e$trip_update$stop_time_update[[1]]
                    tibble(vehicle_id = e$trip_update$vehicle$id,
                           trip_id = e$trip_update$trip$trip_id,
                           route_id = e$trip_update$trip$route_id,
                           timestamp = as.POSIXct(e$trip_update$timestamp, origin = "1970-01-01"),
                           stop_sequence = stu$stop_sequence,
                           type = ifelse(stu$has('arrival'), "arrival", "departure"),
                           time = as.POSIXct(ifelse(stu$has('arrival'), stu$arrival$time, stu$departure$time), 
                                             origin = "1970-01-01")
                    )
                })
            )
        })
    arrivaldata <- do.call(bind_rows, arrivaldata) %>%
        group_by(vehicle_id, route_id, trip_id, stop_sequence, type) %>%
        summarize(time = tmax(timestamp, time))
    save(arrivaldata, file = "simulations/arrivaldata.rda")
}


ids <- tapply(arrivaldata$trip_id, arrivaldata$route_id, function(x) sort(unique(x)))
arrivaldata <- arrivaldata %>% group_by(trip_id)

## get some info from the database
con <- dbConnect(SQLite(), "fulldata.db")
routes <- con %>% tbl("routes") %>% select(route_id, route_short_name, route_long_name) %>% collect
trips <- con %>% tbl("stop_times") %>% filter(stop_sequence == 1) %>% select(trip_id, departure_time) %>% collect
dbDisconnect(con)

get_stop_times <- function(tid) {
    con <- dbConnect(SQLite(), "fulldata.db")
    on.exit(dbDisconnect(con))
    con %>% tbl("stop_times") %>% filter(trip_id == tid) %>%
        select(stop_sequence, arrival_time) %>% arrange(stop_sequence) %>% collect
}

routes <- routes[routes$route_id %in% names(ids), ] %>%
    arrange(route_short_name, route_long_name)
trips <- trips[trips$trip_id %in% arrivaldata$trip_id, ] %>% 
    arrange(departure_time)
rl <- as.list(routes$route_id)
names(rl) <- paste(routes$route_short_name, routes$route_long_name)

## Clean them at the start of the vis
# system("rm -f simulations/*/etas/*.rda")

astime <- function(x) as.POSIXct(x, origin = "1970-01-01")

## Load simulations
whatOrder <- c("number of vehicles", "loading vehicle positions", 
               "updating vehicle information", "updating vehicle states", 
               "predicting ETAs", "writing ETAs to protobuf feed")
view <- function(sim) {
    config <- jsonlite::read_json(file.path("simulations", sim, "config.json"))
    timings <- read_csv(file.path("simulations", sim, "timings.csv"))
    timings <- timings %>% bind_rows(timings %>% group_by(iteration) %>% 
            summarize(timestamp = mean(timestamp), nvehicles = mean(nvehicles)) %>%
            mutate(what = "number of vehicles",
                   cpu = nvehicles, wall = nvehicles)
            ) %>%
        mutate(what = fct_relevel(what, whatOrder),
               timestamp = as.POSIXct(timestamp, origin = "1970-01-01"))

    ggplot(timings, aes(timestamp, cpu)) +
        geom_line() +
        facet_grid(what~., scales = "free_y")
}
eta <- function(sim, ta, ts, sched) {
    config <- jsonlite::read_json(file.path("simulations", sim, "config.json"))
    
    etatime <- as.POSIXct(as.numeric(ts), origin = "1970-01-01")
    etadata <- loadsim(sim, as.integer(ts))
    # cat(" 3 ------\n")
    # print(etadata)
    # cat(" 4 ------\n")
    # print(unique(etadata$trip_id))
    # cat(" 5 ------\n")
    # print(ta$trip_id[1])
    etadata <- etadata %>% filter(trip_id == ta$trip_id[1])
    # cat(" 6 ------\n")
    # print(etadata %>% select(timestamp, stop_sequence, time))
    
    ## calculate the actual arrival/departure time at the most recently visited stop
    delay <- 0
    Slast <- ta %>% filter(time <= astime(ts)) %>% tail(1)
    if (nrow(Slast) == 1)
        delay <- Slast$time - sched$arrival_time[Slast$stop_sequence]

    if (nrow(Slast) == 1)
        sched2 <- sched %>% filter(stop_sequence > Slast$stop_sequence)
    else sched2 <- sched
    p <- ggplot(ta, aes(time, stop_sequence)) +
        geom_point(aes(colour = type)) +
        geom_point(aes(x = arrival_time), data = sched, color = "magenta", pch = 4) +
        geom_point(aes(x = arrival_time + delay), data = sched2, color = "magenta", pch = 2) +
        geom_vline(aes(xintercept = etatime), data = NULL, col = "red", lty = 3) + 
        ggtitle(sprintf("ETAs at %s", format(etatime, "%H:%M:%S"))) +
        xlim(min(ta$time, sched$arrival_time), max(ta$time, sched$arrival_time))
    # print(etadata %>% select(stop_sequence, time, q5, q100))
    if (all(c("q5", "q95") %in% names(etadata))) {
        p <- p + geom_segment(aes(x = astime(q5), xend = pmin(max(ta$time, sched$arrival_time), astime(q95)), yend = stop_sequence), 
            data = etadata %>% filter(q5 > 0 & q95 > 0))
    }
    if ("q50" %in% names(etadata)) {
        p <- p + geom_point(aes(x = astime(q50)), data = etadata %>% filter(q50 > 0), color = "black")
    } else {
        p <- p + geom_point(data = etadata, color = "black", pch = 3)
    }
    p
}
vehicle <- function(vps, ts, prop) {
    obstime <- as.POSIXct(as.numeric(ts), origin = "1970-01-01")
    vps <- vps %>% 
        mutate(dist_between = geosphere::distGeo(cbind(obs_lon, obs_lat), cbind(model_lon, model_lat)),
               ts = as.POSIXct(timestamp, origin = "1970-01-01"))
    p <- vps %>% filter(timestamp == as.integer(obstime))
    con <- dbConnect(SQLite(), "fulldata.db")
    sid <- con %>% tbl("trips") %>% filter(trip_id == p$trip_id[1]) %>% select(shape_id) %>% collect
    shape <- con %>% tbl("shapes") %>% filter(shape_id == sid$shape_id[1]) %>%
        select(shape_pt_lon, shape_pt_lat) %>% arrange(shape_pt_sequence) %>% collect
    dbDisconnect(con)
    vpx <- vps %>% group_by(timestamp) %>%
        summarize(obs_lon = first(obs_lon), obs_lat = first(obs_lat), trip_id = first(trip_id))
    prop <- prop %>% mutate(ts = as.POSIXct(timestamp, origin = "1970-01-01"))
    Dmax <- prop %>% filter(trip_id == p$trip_id) %>% pluck("distance") %>% max
    p1 <- ggplot(p, aes(obs_lon, obs_lat)) +
        coord_fixed(1.6) +
        theme(legend.position = 'none') +
        geom_path(aes(shape_pt_lon, shape_pt_lat), data = shape, colour = "cyan") +
        ## proposals
        geom_point(aes(model_lon, model_lat), pch = 19, col = 'gray', 
            data = prop %>% filter(timestamp == as.integer(obstime))) +
        ## particles
        geom_point(aes(model_lon, model_lat), pch = 4) +
        ## observation
        geom_point(aes(colour = trip_id, alpha = timestamp == as.integer(obstime)), 
            data = vpx)
    p2 <- ggplot(vps %>% filter(trip_id == p$trip_id), aes(distance, speed/1000*60*60)) +
        geom_point(pch = 19, col = 'orangered', 
            data = prop %>% filter(timestamp == as.integer(obstime))) +
        geom_point(col = 'gray', alpha = 0.1, pch = 4) +
        geom_point(data = p, alpha = 1) +
        theme(legend.position = 'none') + ylim(0, 100) + ylab("Speed (km/h)") +
        xlab("Distance Traveled (m)") + xlim(0, Dmax)
    p3 <- ggplot(vps %>% filter(trip_id == p$trip_id), aes(ts, distance)) +
        geom_point(pch = 19, col = 'orangered', 
            data = prop %>% filter(timestamp == as.integer(obstime))) +
        geom_point(col = 'gray', alpha = 0.1, pch = 4) +
        geom_point(data = p, alpha = 1) +
        theme(legend.position = 'none') +
        xlab("Time") +
        ylab("Distance Traveled (m)") + ylim(0, Dmax)
    p4 <- ggplot(vps %>% filter(trip_id == p$trip_id) %>% group_by(timestamp) %>%
                summarize(x = mean(distance), xdot = mean(speed)) %>%
                ungroup() %>% mutate(avg_speed = c(0, diff(x) / diff(timestamp)))) +
        geom_point(aes(x, avg_speed/1000*60*60, colour = timestamp == as.integer(obstime))) +
        ylim(0, 100) + ylab("Mean Interpolated Speed (km/h)") +
        xlab("Distance Traveled (m)") + xlim(0, Dmax) +
        theme(legend.position = 'none')
    p5 <- ggplot(vps %>% filter(trip_id == p$trip_id), aes(ts, dist_between)) +
        geom_point(aes(colour = timestamp == as.integer(obstime))) +
        theme(legend.position = 'none') +
        xlab("Time") +
        ylab("Distance between\nparticle and observation (m)")
    p6 <- ggplot(p, aes(distance, ll)) + 
        geom_point(col = 'orangered', data = prop %>% filter(timestamp == as.integer(obstime))) +
        geom_point() +
        theme(legend.position = 'none') +
        xlab("Distance Traveled (m)") + xlim(0, Dmax) +
        ylab("Log Likelihood")

    gridExtra::grid.arrange(p1, p3, p5, p2, p4, p6, 
        layout_matrix = cbind(c(1, 1, 1, 1, 6), c(2:6)))
}

modeleval <- function(data, n, gps) {
    data <- data %>% mutate(ts = astime(ts))

    p <- ggplot(data, aes(ts)) + xlab("Time")

    p_neff <- p + geom_point(aes(y = Neff, color = factor(resample, levels = 0:1))) +
        theme(legend.position = "top") + ylim(0, n) + ylab("Effective Sample Size")
    p_prior <- p + 
        geom_hline(aes(yintercept = gps^2), lty = 3) + 
        geom_path(aes(y = pmin(10 * gps^2, dist_to_path)), lty = 2, col = 'darkblue') + 
        geom_path(aes(y = pmin(10 * gps^2, prior_mse %>% sqrt))) + 
        geom_path(aes(y = pmin(10 * gps^2, posterior_mse %>% sqrt)), color = 'orangered') +
        ylim(0, 10 * gps^2) + ylab("MSE (m)") +
        facet_wrap(~trip_id, scales = "free_x")

    gridExtra::grid.arrange(
        p_neff, p_prior,
        layout_matrix = cbind(1:2), heights = c(1, 3))
}

library(shiny)
ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            selectInput("simnum", label="Simulation", choices = simnames),
            verbatimTextOutput("simconfig"),
            hr(),
            selectInput("plottype", label="Graph", choices = c("Timings", "Vehicles", "ETAs", "Model Evaluation")),
            conditionalPanel(condition = "input.plottype == 'ETAs'",
                selectInput("routeid", label="Route ID", choices = rl),
                selectInput("tripid", label="Trip ID", ""),
                sliderInput("predtime", label="Prediction Time", min = 1, max = 1, step=1, value = 1)
            ),
            conditionalPanel(condition = "input.plottype == 'Vehicles' || input.plottype == 'Model Evaluation'",
                selectInput("vehicleid", label="Vehicle ID", choices = "")
            ),
            conditionalPanel(condition = "input.plottype == 'Vehicles'",
                sliderInput("obstime", label="Observation Time", min = 1, max = 1, step=1, value = 1)
            )
        ),
        mainPanel(
            plotOutput(outputId = "graph", width="800", height="900")
        )
    )
)
server <- function(input, output, session) {
    # refresh <- reactiveTimer(2000)
    pbapply::pboptions(type = "timer")
    rv <- reactiveValues()
    rv$ta <- NULL
    rv$etatimes <- NULL
    rv$tatimes <- NULL
    rv$schedule <- NULL
    rv$config <- NULL
    rv$vehicledata <- NULL
    rv$hdir <- NULL
    observeEvent(input$simnum, {
        ## config info
        conf <- readLines(file.path("simulations", input$simnum, "config.json"))
        output$simconfig <- renderPrint(jsonlite::prettify(conf))

        rv$hdir <- file.path("simulations", input$simnum, "history")
        rv$config <- jsonlite::fromJSON(conf)
        if (!is.null(rv$config$simulation_history)) 
            rv$hdir <- rv$config$simulation_history
        if (dir.exists(rv$hdir)) {
            vids <- gsub("vehicle_|\\.csv", "", list.files(rv$hdir))
            vids <- vids[!grepl("_", vids)]
            updateSelectInput(session, "vehicleid", choices = vids)
        }

        ## sim stuff
        fts <- list.files(file.path("simulations", input$simnum, "etas"), pattern = "*.pb")
        if (length(fts) == 0) {
            rv$etatimes <- NULL
            rv$tatimes <- NULL
        } else {
            # rv$etatimes <- all_sims(input$simnum) %>% filter(!is.na(time))
            # rv$tatimes <- rv$etatimes %>% filter(trip_id == )
            rv$etatimes <- as.POSIXct(as.numeric(gsub("etas_|\\.pb", "", fts)), origin = "1970-01-01")
            rv$tatimes <- rv$etatimes[rv$etatimes >= min(rv$ta$time) &
                                      rv$etatimes <= max(rv$ta$time)]
        }
    })
    observeEvent(input$routeid, {
        ti <- trips[trips$trip_id %in% ids[[input$routeid]], ]
        tl <- as.list(ti$trip_id)
        names(tl) <- ti$departure_time
        updateSelectInput(session, "tripid", choices = tl)
    })
    observeEvent(input$tripid, {
        if (input$tripid != "") {
            rv$ta <- arrivaldata %>% filter(trip_id == input$tripid)
            day <- format(rv$ta$time[1], "%Y-%m-%d")
            
            tstart <- as.POSIXct(paste(day, trips[trips$trip_id == rv$ta$trip_id[1], "departure_time"]))
            rv$ta <- rv$ta %>% filter(stop_sequence > 1 | time < tstart + 60*30)
            rv$schedule <- get_stop_times(input$tripid) %>% 
                mutate(arrival_time = as.POSIXct(paste(day, arrival_time)))
            if (!is.null(rv$etatimes)) {
                # rv$tatimes <- rv$etatimes %>% filter(trip_id == input$tripid) %>%
                #     pluck("timestamp") %>% unique() %>% sort()
                rv$tatimes <- rv$etatimes[rv$etatimes >= min(rv$ta$time) - 60*30 &
                                          rv$etatimes <= max(rv$ta$time) + 60*60]
            } else {
                rv$tatimes <- NULL
                rv$tatimes <- NULL
            }
            updateSliderInput(session, "predtime", max = length(rv$tatimes))
        }
    })
    observeEvent(input$vehicleid, {
        if (input$vehicleid == "") return()
        if (input$plottype == "Vehicles") {
            rv$vehicledata <- read_csv(sprintf("%s/vehicle_%s.csv", rv$hdir, input$vehicleid),
                col_names = c('timestamp', 'trip_id', 'obs_lat', 'obs_lon', 'distance', 'speed', 'acceleration', 'll', 'model_lat', 'model_lon'))
            rv$vehicleprop <- read_csv(sprintf("%s/vehicle_%s_proposals.csv", rv$hdir, input$vehicleid),
                col_names = c('timestamp', 'trip_id', 'distance', 'speed', 'acceleration', 'll', 'model_lat', 'model_lon'))
            rv$vehicletimes <- unique(rv$vehicledata$timestamp)
            updateSliderInput(session, "obstime", max = length(rv$vehicletimes))
        } else {
            rv$vehicledata <- read_csv(sprintf("simulations/%s/modeleval/vehicle_%s.csv", input$simnum, input$vehicleid),
                col_names = c("vehicle_id", "trip_id", "ts", "prior_mse", "posterior_mse", "dist_to_path", "Neff", "resample"))
            # print(rv$vehicledata)
        }
    })
    observeEvent(input$plottype, {
        if (input$plottype == "Vehicles" && 
            dir.exists(rv$hdir)) {
            vids <- gsub("vehicle_|\\.csv", "", list.files(rv$hdir))
            vids <- vids[!grepl("_", vids)]
            updateSelectInput(session, "vehicleid", choices = vids)
        }
        if (input$plottype == "Model Evaluation" && 
            dir.exists(file.path("simulations", input$simnum, "modeleval"))) {
            vids <- gsub("vehicle_|\\.csv", "", list.files(file.path("simulations", input$simnum, "modeleval")))
            updateSelectInput(session, "vehicleid", choices = vids)
        }    
    })
    output$graph <- renderPlot({
        # refresh()
        switch(input$plottype, 
            "Timings" = view(input$simnum),
            "ETAs" = eta(input$simnum, rv$ta, rv$tatimes[input$predtime], rv$schedule),
            "Vehicles" = vehicle(rv$vehicledata, rv$vehicletimes[input$obstime], rv$vehicleprop),
            "Model Evaluation" = modeleval(rv$vehicledata, rv$config$n_particles, rv$config$gps_error)
        )
    })
}
shinyApp(ui = ui, server = server, options = list(port = 9999))
