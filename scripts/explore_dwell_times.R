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

# library(biglm)

# i <- 1
# size <- 1000
# fit <- biglm(dwell ~ stop_id, data = dt[i:(i+size-1),])
# i <- i+size
# cat("\n")
# while (i < nrow(dt)) {
#     cat(sprintf("\r%i / %i ", i, nrow(dt)))
#     fit <- update(fit, dt[i:(i+size-1),])
#     i <- i+size
# }
# cat("\r done                  \n")



### Look at travel times

trip_times_to_travel_times <- function(x) {
    ## sort & remove duplicates
    x <- x %>% distinct() %>% arrange(stop_sequence, type) %>%
        mutate(time = as.POSIXct(time, origin = "1970-01-01"))
    ## all should be increasing
    if (nrow(x) < 3 || any(diff(x$time) < 0)) 
        return(x %>% filter(FALSE) %>% 
            mutate(segment_index = stop_sequence, 
                time_from = time, time_to = time, travel_time = numeric()) %>%
            select(route_id, trip_id, segment_index, time_from, time_to, travel_time))
    ## rearrange into from-to
    xFrom <- x %>% 
        filter(stop_sequence < max(stop_sequence) & type == "departure") %>%
        mutate(segment_index = stop_sequence)
    xTo <- x %>% filter(stop_sequence > 1 & type == "arrival") %>%
        mutate(segment_index = stop_sequence - 1)
    xFrom %>% select(route_id, trip_id, segment_index, time) %>%
        inner_join(
            xTo %>% select(segment_index, time), 
            by = "segment_index", 
            suffix = c("_from", "_to")
        ) %>% 
        mutate(
            time_from = time_from,
            time_to = time_to,
            travel_time = as.numeric(time_to - time_from)
        )
}
# x <- tu %>% filter(trip_id == tu$trip_id[1])
# x <- x %>% filter(vehicle_id == x$vehicle_id %>% table %>% sort %>% tail(1) %>% names)

tu0 <- tu %>% distinct() %>%
    group_by(trip_id, vehicle_id) %>%
    do((.) %>% trip_times_to_travel_times())

library(RSQLite)
library(dbplyr)
rids <- unique(tu0$route_id)
con <- dbConnect(SQLite(), "fulldata.db")
routei <- con %>% tbl("routes") %>%
    filter(route_id %in% rids) %>%
    select(route_id, route_short_name, route_long_name) %>%
    collect()
dbDisconnect(con)
rnums <- structure(routei$route_short_name, .Names = routei$route_id)

tu1 <- tu0 %>% mutate(route_number = rnums[route_id])

    #     filter(between(
    #         timestamp, 
    #         median(timestamp) - 60 * 60, 
    #         median(timestamp) + 60 * 60
    #     ))
    # ) %>%
    # ungroup() %>%
    # group_by(vehicle_id, route_id, trip_id, timestamp) %>%
    # arrange(timestamp) %>%
    # summarize(
    #     stop_sequence = first(stop_sequence),
    #     type = first(type),
    #     time = first(time)
    # )

std <- function(x) {
    if (length(x) == 1) return(0)
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

rids <- names(sort(table(tu$route_id), TRUE))
RID <- rids[10]; routei %>% filter(route_id == RID)

tu1r <- tu1 %>% 
    filter(
        route_id == RID &
        between(travel_time, 20, 60*30)
    ) %>%
    filter(
        segment_index < max((.)$segment_index)
    ) %>%
    group_by(segment_index) %>%
    do((.) %>% mutate(
        travel_time.std = std(travel_time)
    )) %>%
    ungroup() %>%
    arrange(segment_index, time_from)
 

ggplot(
    tu1r, 
    aes(time_from, travel_time/60)
) + geom_point() + facet_wrap(~segment_index)

ggplot(
    tu1r, 
    aes(time_from, travel_time.std/60)
) + geom_point() + facet_wrap(~segment_index)

## num -> [30min intervals]
times <- pretty(tu1r$time_from, 20 * 2)
smry <- tu1r %>%
    mutate(
        period = cut((.)$time_from, breaks = times)
    ) %>%
    group_by(segment_index, period) %>%
    summarize(
        tt = mean(travel_time, na.rm = TRUE),
        sd.tt = sd(travel_time, na.rm = TRUE),
        n = n()
    ) %>%
    mutate(
        time = as.POSIXct(period)
    )

ggplot(smry, aes(time, tt/60)) + 
    geom_point() +
    facet_wrap(~segment_index)






## now a combination of [trip_id, shape_id, segment_index] can get [segment_id]
con <- dbConnect(SQLite(), "fulldata.db")
tids <- unique(tu1$trip_id)
tsegs <- con %>% tbl("trips") %>% 
    filter(trip_id %in% tids) %>%
    left_join(con %>% tbl("shape_segments")) %>%
    left_join(con %>% tbl("road_segments") %>% select(road_segment_id, length)) %>%
    select(trip_id, road_segment_id, shape_road_sequence, length) %>%
    collect()
segids <- unique(tsegs$road_segment_id)
segs <- con %>% tbl("road_segments") %>%
    filter(road_segment_id %in% segids) %>%
    left_join(con %>% tbl("intersections"), by = c("int_from" = "intersection_id")) %>%
    left_join(con %>% tbl("intersections"), by = c("int_to" = "intersection_id"), suffix = c("", "_to")) %>%
    collect()
dbDisconnect(con)

sx <- tu1 %>% left_join(tsegs, by = c("trip_id", "segment_index" = "shape_road_sequence")) %>%
    filter(travel_time > 20 & travel_time < 60*30) %>%
    arrange(road_segment_id, time_from) %>%
    mutate(speed = length / travel_time) %>%
    filter(speed < 30)

seg20 <- sx$road_segment_id %>% table %>% sort(TRUE) %>% names
sx20 <- sx %>% filter(road_segment_id %in% seg20[61:80])
# ggplot(sx20, aes(time_from, speed / 1000 * 60 * 60)) + 
ggplot(sx20, aes(time_from, travel_time / 60)) + 
    geom_point(aes(colour = route_id)) +
    facet_wrap(~road_segment_id, scales = "free_y") +
    xlab("Time") + 
    # ylab("Speed (km/h)") +
    ylab("Travel time (min)") +
    scale_colour_discrete(guide = FALSE) + 
    scale_x_datetime(label = function(x) format(x, '%H:%M'))
    # scale_y_log10()
    # geom_smooth(method = lm, formula = y ~ splines::bs(x, 4))


### I want the model to work on both mu and sigma
# i.e., mu(t) and sigma(t)

r1 <- 







# tu0s <- tu0 %>% ungroup() %>%
#     group_by(vehicle_id, route_id, trip_id, stop_sequence) %>%
#     summarize(
#         arr = min(time),
#         dep = max(time),
#         dwell = dep - arr
#     ) %>%
#     do((.) %>%
#         arrange(stop_sequence) %>%
#         mutate(travel_time = c(0, .$arr[-1] - .$dep[-n()]))
#     )
# ## arr time at 'next' stop - dep time at 'this' stop

# ggplot(tu0s, aes(trip_id, travel_time)) +
#     geom_point()

# trs <- tu0s %>% group_by(vehicle_id, route_id, trip_id) %>%
#     summarize(
#         time = min(arr),
#         trip_time = sum(travel_time),
#         n_stops = n(),
#         min_stop = min(stop_sequence),
#         max_stop = max(stop_sequence)
#     ) 

# ## only those with full obs [trip_id -> stop_count]
# library(RSQLite)
# library(dbplyr)
# tids <- unique(trs$trip_id)
# con <- dbConnect(SQLite(), "fulldata.db")
# stopn <- con %>% tbl("stop_times") %>% 
#     filter(trip_id %in% tids) %>%
#     group_by(trip_id) %>%
#     summarize(n = n()) %>%
#     collect()
# dbDisconnect(con)

# stopn <- structure(stopn$n, .Names = stopn$trip_id)
# trs2 <- trs %>% filter(min_stop == 1 & max_stop == stopn[trip_id])

# ggplot(trs, aes(time, trip_time/60)) + geom_point()


# ## can we map each observed segment to a road_segment ???
