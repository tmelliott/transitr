library(tidyverse)
library(RProtoBuf)

curd <- setwd("src/vendor/protobuf")
readProtoFiles("gtfs-realtime.proto")
setwd(curd)

## load trip updates
if (!file.exists("tripupdates.rda")) {
    tufiles <- list.files(file.path("simulations", "archive"), 
        pattern = "trip_updates*", full.name = TRUE)
    tu <- pbapply::pblapply(tufiles, function(f) {
        feed <- read(transit_realtime.FeedMessage, f)
        if (length(feed$entity) == 0) return(NULL)
        lapply(feed$entity, function(e) {
            stus <- e$trip_update$stop_time_update
            if (length(stus) == 0) return(NULL)
            xdf <- tibble(
                vehicle_id = e$trip_update$vehicle$id,
                trip_id = e$trip_update$trip$trip_id,
                route_id = e$trip_update$trip$route_id,
                timestamp = as.POSIXct(e$trip_update$timestamp, origin = "1970-01-01"),
                stop_sequence = sapply(stus, function(stu) 
                    if (stu$has('stop_sequence')) stu$stop_sequence else NA
                ),
                type = sapply(stus, function(stu) ifelse(stu$has('arrival'), 'arrival', 'departure')),
                time = sapply(stus, function(stu) if (stu$has('arrival')) stu$arrival$time else stu$departure$time)
            )
        }) %>% bind_rows
    }) %>% bind_rows
    save(tu, file = "tripupdates.rda")
} else {
    load("tripupdates.rda")
}

## load vehilce positions
if (!file.exists("vehiclepositions.rda")) {
    vpfiles <- list.files(file.path("simulations", "archive"), 
        pattern = "vehicle_locations*", full.name = TRUE)
    vp <- pbapply::pblapply(vpfiles, function(f) {
        feed <- read(transit_realtime.FeedMessage, f)
        if (length(feed$entity) == 0) return(NULL)
        lapply(feed$entity, function(e) {
            if (!e$has('vehicle')) return(NULL)
            if (!e$vehicle$has('position')) return(NULL)
            xdf <- tibble(
                vehicle_id = e$vehicle$vehicle$id,
                trip_id = e$vehicle$trip$trip_id,
                route_id = e$vehicle$trip$route_id,
                timestamp = as.POSIXct(e$vehicle$timestamp, origin = "1970-01-01"),
                position_latitude = e$vehicle$position$latitude,
                position_longitude = e$vehicle$position$longitude
            )
        }) %>% bind_rows
    }) %>% bind_rows
    save(vp, file = "vehiclepositions.rda")
} else {
    load("vehiclepositions.rda")
}

dt <- tu %>%
    filter(type == "arrival") %>%
    mutate(arrival = time) %>%
    select(vehicle_id, trip_id, timestamp, stop_sequence, arrival) %>%
    distinct() %>% 
    full_join(
        tu %>% filter(type == "departure") %>% mutate(departure = time) %>%
            select(vehicle_id, trip_id, stop_sequence, departure),
        # by = c("vehicle_id", "trip_id", "stop_sequence"),
        suffix = c(".arrival", ".departure")
    ) %>%
    arrange(vehicle_id, trip_id, stop_sequence) %>%
    mutate(dwell = departure - arrival) %>% filter(!is.na(dwell)) %>%
    ## only include dwell times [0, 1min]
    filter(dwell >= 0 & dwell <= 1*60)

hp <- ggplot(dt, aes(dwell)) + 
    geom_histogram(binwidth = 1) + 
    xlab("Dwell time (sec)") + ylab("Number of buses") 

hp + geom_vline(aes(xintercept = x), 
    data = tibble(x = c(9, 19, 27, 36, 45)), color = "orangered")

hp + geom_vline(aes(xintercept = x), 
    data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) 


## how about the rate of vehicle updates?
dv <- vp %>% group_by(vehicle_id) %>% distinct() %>%
    do((.) %>% arrange(timestamp) %>% mutate(dt = c(-1, diff(timestamp)))) %>%
    filter(dt >= 0 & dt <= 1*60)

dp <- ggplot(dv, aes(dt)) +
    geom_histogram(binwidth = 1) +
    xlab("Time between updates (sec)") + ylab("Number of updates")

dp + geom_vline(aes(xintercept = x), 
    data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) 

dt2 <- tu %>% group_by(vehicle_id) %>% distinct() %>%
    do((.) %>% arrange(timestamp) %>% mutate(dt = c(-1, diff(timestamp)))) %>%
    filter(dt >= 0 & dt <= 1*60)

hp2 <- ggplot(dt2, aes(dt)) +
    geom_histogram(binwidth = 1) +
    xlab("Time between updates (sec)") + ylab("Number of updates")

hp2 + geom_vline(aes(xintercept = x), 
    data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) 


## join
dd <- dv %>% full_join(dt, by = c("vehicle_id", "timestamp")) %>%
    distinct() %>%
    do((.) %>% arrange(timestamp) %>% mutate(dt = c(-1, diff(timestamp)))) %>%
    filter(dt >= 0 & dt <= 1*60)

dp2 <- ggplot(dd %>% filter(dt < 1*60), aes(dt)) +
    geom_histogram(binwidth = 1) +
    xlab("Time between updates (sec)") + ylab("Number of updates")

dp2 + geom_vline(aes(xintercept = x), 
    data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) 

gridExtra::grid.arrange(
    hp + geom_vline(aes(xintercept = x), 
        data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) +
        ggtitle("Dwell times"),
    hp2 + geom_vline(aes(xintercept = x), 
        data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) +
        ggtitle("Trip updates only"),
    dp + geom_vline(aes(xintercept = x), 
        data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) +
        ggtitle("Vehicle positions only"),
    dp2 + geom_vline(aes(xintercept = x), 
        data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) +
        ggtitle("Combined"),
    nrow = 4)



## find the stop
library(RSQLite)
library(dbplyr)
tids <- unique(dt$trip_id)
sts <- dbConnect(SQLite(), "fulldata.db") %>% tbl("stop_times") %>% 
    select(trip_id, stop_sequence, stop_id) %>%
    filter(trip_id %in% tids) %>%
    collect

dt <- dt %>% left_join(sts) %>%
    mutate(
        stop_id = as.factor(stop_id),
        stopid = as.numeric(stop_id)
    )

library(biglm)

i <- 1
size <- 1000
fit <- biglm(dwell ~ stop_id, data = dt[i:(i+size-1),])
i <- i+size
cat("\n")
while (i < nrow(dt)) {
    cat(sprintf("\r%i / %i ", i, nrow(dt)))
    fit <- update(fit, dt[i:(i+size-1),])
    i <- i+size
}
cat("\r done                  \n")
