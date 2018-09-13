library(tidyverse)
library(RProtoBuf)

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
loadsim <- function(sim) {
    if (file.exists(file.path("simulations", sim, "etadata.rda"))) {
        load(file.path("simulations", sim, "etadata.rda"))
    } else {
        etas <- do.call(bind_rows, 
            pbapply::pblapply(list.files(file.path("simulations", sim, "etas"), full.names = TRUE),
                function(f) {
                    ent <- read(transit_realtime.FeedMessage, f)$entity
                    do.call(bind_rows, 
                        lapply(ent, function(e) {
                            lapply(e$trip_update$stop_time_update, function(stu) {
                                eta <- stu$getExtension(transit_network.eta)
                                tibble(vehicle_id = e$trip_update$vehicle$id,
                                       trip_id = e$trip_update$trip$trip_id,
                                       route_id = e$trip_update$trip$route_id,
                                       timestamp = as.POSIXct(as.numeric(gsub("^.+_|\\.pb$", "", f)), origin = "1970-01-01"),
                                       stop_sequence = stu$stop_sequence,
                                       time = as.POSIXct(ifelse(eta$estimate == 0, NA, eta$estimate), origin = "1970-01-01")
                                )
                            })
                        })
                    )
                })
        )
        save(etas, file = file.path("simulations", sim, "etadata.rda"))
    }
    etas
}

ids <- tapply(arrivaldata$trip_id, arrivaldata$route_id, function(x) sort(unique(x)))
arrivaldata <- arrivaldata %>% group_by(trip_id)

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
eta <- function(sim, ta, eta) {
    config <- jsonlite::read_json(file.path("simulations", sim, "config.json"))
    
    print(ta)
    print(eta)
    ggplot(ta, aes(time, stop_sequence, colour = type)) +
        geom_point() +
        geom_point(data = eta, color = "red")
}

library(shiny)
ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            selectInput("simnum", label="Simulation", choices = simnames),
            selectInput("plottype", label="Graph", choices = c("Timings", "ETAs")),
            conditionalPanel(condition = "input.plottype == 'ETAs'",
                selectInput("routeid", label="Route ID", choices = names(ids)),
                selectInput("tripid", label="Trip ID", ""),
                sliderInput("predtime", label="Prediction Time", min = 0, max = 1, step=1, value = 0)
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
    rv$etat <- NULL
    rv$etas <- NULL
    observeEvent(input$simnum, {
        rv$etas <- loadsim(input$simnum)
        print(rv$etas)
    })
    
    observeEvent(input$routeid, {
        updateSelectInput(session, "tripid", 
            choices = ids[[input$routeid]])
    })
    observeEvent(input$tripid, {
        rv$ta <- arrivaldata %>% filter(trip_id == input$tripid)
        print(rv$ta)
        rv$etat <- rv$etas %>% filter(trip_id == input$tripid)
        print(rv$etat)
        updateSliderInput(session, "predtime",
            max = length(unique(rv$etat$timestamp)))
    })
    output$graph <- renderPlot({
        # refresh()
        switch(input$plottype, 
            "Timings" = view(input$simnum),
            "ETAs" = eta(input$simnum, rv$ta, rv$etat)
        )
    })
}
shinyApp(ui = ui, server = server, options = list(port = 9999))
