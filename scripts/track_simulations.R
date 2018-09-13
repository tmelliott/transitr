library(tidyverse)
library(RProtoBuf)
library(RSQLite)
library(dbplyr)

curd <- setwd("src/vendor/protobuf")
readProtoFiles("gtfs-realtime-ext.proto")
setwd(curd)

simnames <- list.files("simulations", pattern = "sim*", include.dirs = TRUE)
if (file.exists("simulations/arrivaldata.rda")) {
    load("simulations/arrivaldata.rda")
} else {
    tmax <- function(x, t) x[which.max(t)]
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
loadsim <- function(sim, time) {
    pb <- file.path("simulations", sim, "etas", sprintf("etas_%s.pb", time))
    rda <- file.path("simulations", sim, "etas", sprintf("etas_%s.rda", time))
    if (!file.exists(pb)) {
        stop("That file doesn't exist.")
    }
    if (file.exists(rda)) {
        load(rda)
    } else {
        ent <- read(transit_realtime.FeedMessage, pb)$entity
        etas <- do.call(bind_rows, 
            pbapply::pblapply(ent, function(e) {
                lapply(e$trip_update$stop_time_update, function(stu) {
                    eta <- stu$getExtension(transit_network.eta)
                    tibble(vehicle_id = e$trip_update$vehicle$id,
                           trip_id = e$trip_update$trip$trip_id,
                           route_id = e$trip_update$trip$route_id,
                           timestamp = as.POSIXct(as.integer(time), origin = "1970-01-01"),
                           stop_sequence = if (stu$has('stop_sequence')) stu$stop_sequence else NA,
                           time = if (eta$has('estimate') && eta$estimate > 0) eta$estimate else NA
                    )
                })
            })
        ) %>%
            mutate(time = as.POSIXct(time, origin = "1970-01-01")) %>%
            group_by(trip_id)
        save(etas, file = rda)
    }
    etas
}

ids <- tapply(arrivaldata$trip_id, arrivaldata$route_id, function(x) sort(unique(x)))
arrivaldata <- arrivaldata %>% group_by(trip_id)

## get some info from the database
con <- dbConnect(SQLite(), "fulldata.db")
routes <- con %>% tbl("routes") %>% select(route_id, route_short_name, route_long_name) %>% collect
trips <- con %>% tbl("stop_times") %>% filter(stop_sequence == 1) %>% select(trip_id, departure_time) %>% collect
dbDisconnect(con)

routes <- routes[routes$route_id %in% names(ids), ] %>%
    arrange(route_short_name, route_long_name)
trips <- trips[trips$trip_id %in% arrivaldata$trip_id, ] %>% 
    arrange(departure_time)
rl <- as.list(routes$route_id)
names(rl) <- paste(routes$route_short_name, routes$route_long_name)

## Clean them at the start of the vis
# system("rm -f simulations/*/etas/*.rda")

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
eta <- function(sim, ta, ts) {
    config <- jsonlite::read_json(file.path("simulations", sim, "config.json"))
    
    etatime <- as.POSIXct(as.numeric(ts), origin = "1970-01-01")
    etadata <- loadsim(sim, as.integer(ts)) %>%
        filter(trip_id == ta$trip_id[1])
    ggplot(ta, aes(time, stop_sequence)) +
        geom_point(aes(colour = type)) +
        geom_vline(aes(xintercept = etatime), data = NULL, col = "red", lty = 3) + 
        ggtitle(sprintf("ETAs at %s", format(etatime, "%H:%M:%S"))) +
        geom_point(data = etadata, color = "black", pch = 3)
}

library(shiny)
ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            selectInput("simnum", label="Simulation", choices = simnames),
            selectInput("plottype", label="Graph", choices = c("Timings", "ETAs")),
            conditionalPanel(condition = "input.plottype == 'ETAs'",
                selectInput("routeid", label="Route ID", choices = rl),
                selectInput("tripid", label="Trip ID", ""),
                sliderInput("predtime", label="Prediction Time", min = 1, max = 1, step=1, value = 1)
            )
        ),
        mainPanel(
            plotOutput(outputId = "graph", width="800", height="900")
        )
    )
)
server <- function(input, output, session) {
    # refresh <- reactiveTimer(2000)
    rv <- reactiveValues()
    rv$ta <- NULL
    rv$etatimes <- NULL
    rv$tatimes <- NULL
    observeEvent(input$simnum, {
        fts <- list.files(file.path("simulations", input$simnum, "etas"), pattern = "*.pb")
        if (length(fts) == 0) {
            rv$etatimes <- NULL
            rv$tatimes <- NULL
        } else {
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
            if (!is.null(rv$etatimes)) {
                rv$tatimes <- rv$etatimes[rv$etatimes >= min(rv$ta$time) &
                                          rv$etatimes <= max(rv$ta$time)]
            } else {
                rv$tatimes <- NULL
            }
            updateSliderInput(session, "predtime", max = length(rv$tatimes))
        }
    })
    output$graph <- renderPlot({
        # refresh()
        switch(input$plottype, 
            "Timings" = view(input$simnum),
            "ETAs" = eta(input$simnum, rv$ta, rv$tatimes[input$predtime])
        )
    })
}
shinyApp(ui = ui, server = server, options = list(port = 9999))
