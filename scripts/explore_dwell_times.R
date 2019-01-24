library(tidyverse)
library(RProtoBuf)

# curd <- setwd("src/vendor/protobuf")
# readProtoFiles("gtfs-realtime.proto")
# setwd(curd)

## load trip updates
if (!file.exists("tripupdates.rda")) {
    tufiles <- list.files(pattern = "trip_updates_.*\\.rda")
    tripupdates <- NULL
    for (f in tufiles) {
        load(f)
        tripupdates <- tripupdates %>% bind_rows(tu %>% distinct())
        rm(tu)
    }
    save(tripupdates, file = "tripupdates.rda")
} else {
    load("tripupdates.rda")
}

## load vehilce positions
# if (!file.exists("vehiclepositions.rda")) {
#     vpfiles <- list.files(file.path("simulations", "archive"), 
#         pattern = "vehicle_locations*", full.name = TRUE)
#     vp <- pbapply::pblapply(vpfiles, function(f) {
#         feed <- read(transit_realtime.FeedMessage, f)
#         if (length(feed$entity) == 0) return(NULL)
#         lapply(feed$entity, function(e) {
#             if (!e$has('vehicle')) return(NULL)
#             if (!e$vehicle$has('position')) return(NULL)
#             xdf <- tibble(
#                 vehicle_id = e$vehicle$vehicle$id,
#                 trip_id = e$vehicle$trip$trip_id,
#                 route_id = e$vehicle$trip$route_id,
#                 timestamp = as.POSIXct(e$vehicle$timestamp, origin = "1970-01-01"),
#                 position_latitude = e$vehicle$position$latitude,
#                 position_longitude = e$vehicle$position$longitude
#             )
#         }) %>% bind_rows
#     }) %>% bind_rows
#     save(vp, file = "vehiclepositions.rda")
# } else {
#     load("vehiclepositions.rda")
# }

dt <- tripupdates %>%
    filter(type == "arrival") %>%
    mutate(arrival = time) %>%
    select(vehicle_id, trip_id, timestamp, stop_sequence, arrival) %>%
    distinct() %>% 
    full_join(
        tripupdates %>% filter(type == "departure") %>% mutate(departure = time) %>%
            select(vehicle_id, trip_id, stop_sequence, departure),
        # by = c("vehicle_id", "trip_id", "stop_sequence"),
        suffix = c(".arrival", ".departure")
    ) %>%
    arrange(vehicle_id, trip_id, stop_sequence) %>%
    mutate(dwell = departure - arrival) %>% filter(!is.na(dwell)) %>%
    filter(dwell >= 0 & dwell <= 1*60)

hp <- ggplot(dt, aes(dwell)) + 
    geom_histogram(binwidth = 1) + 
    xlab("Dwell time (sec)") + ylab("Number of buses") 

hp + geom_vline(aes(xintercept = x), 
    data = tibble(x = c(9, 19, 27, 36, 45)), color = "orangered")

hp + geom_vline(aes(xintercept = x), 
    data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) 


## how about the rate of vehicle updates?
# dv <- vp %>% group_by(vehicle_id) %>% distinct() %>%
#     do((.) %>% arrange(timestamp) %>% mutate(dt = c(-1, diff(timestamp)))) %>%
#     filter(dt >= 0 & dt <= 1*60)

# dp <- ggplot(dv, aes(dt)) +
#     geom_histogram(binwidth = 1) +
    # xlab("Time between updates (sec)") + ylab("Number of updates")

# dp + geom_vline(aes(xintercept = x), 
#     data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) 

# dt2 <- tu %>% group_by(vehicle_id) %>% distinct() %>%
#     do((.) %>% arrange(timestamp) %>% mutate(dt = c(-1, diff(timestamp)))) %>%
#     filter(dt >= 0 & dt <= 1*60)

# hp2 <- ggplot(dt2, aes(dt)) +
#     geom_histogram(binwidth = 1) +
#     xlab("Time between updates (sec)") + ylab("Number of updates")

# hp2 + geom_vline(aes(xintercept = x), 
#     data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) 


## join
# dd <- dv %>% full_join(dt, by = c("vehicle_id", "timestamp")) %>%
#     distinct() %>%
#     do((.) %>% arrange(timestamp) %>% mutate(dt = c(-1, diff(timestamp)))) %>%
#     filter(dt >= 0 & dt <= 1*60)

# dp2 <- ggplot(dd %>% filter(dt < 1*60), aes(dt)) +
#     geom_histogram(binwidth = 1) +
#     xlab("Time between updates (sec)") + ylab("Number of updates")

# dp2 + geom_vline(aes(xintercept = x), 
#     data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) 

# gridExtra::grid.arrange(
#     hp + geom_vline(aes(xintercept = x), 
#         data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) +
#         ggtitle("Dwell times"),
#     hp2 + geom_vline(aes(xintercept = x), 
#         data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) +
#         ggtitle("Trip updates only"),
#     dp + geom_vline(aes(xintercept = x), 
#         data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) +
#         ggtitle("Vehicle positions only"),
#     dp2 + geom_vline(aes(xintercept = x), 
#         data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) +
#         ggtitle("Combined"),
#     nrow = 4)



# ## find the stop
# library(RSQLite)
# library(dbplyr)
# tids <- unique(dt$trip_id)
# sts <- dbConnect(SQLite(), "fulldata.db") %>% tbl("stop_times") %>% 
#     select(trip_id, stop_sequence, stop_id) %>%
#     filter(trip_id %in% tids) %>%
#     collect

# dt <- dt %>% left_join(sts) %>%
#     mutate(
#         stop_id = as.factor(stop_id),
#         stopid = as.numeric(stop_id)
#     )

# # library(biglm)

# # i <- 1
# # size <- 1000
# # fit <- biglm(dwell ~ stop_id, data = dt[i:(i+size-1),])
# # i <- i+size
# # cat("\n")
# # while (i < nrow(dt)) {
# #     cat(sprintf("\r%i / %i ", i, nrow(dt)))
# #     fit <- update(fit, dt[i:(i+size-1),])
# #     i <- i+size
# # }
# # cat("\r done                  \n")



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

tu <- tripupdates
rm(tripupdates)
if (file.exists("tu0.rda")) {
    load("tu0.rda")
} else {
    tu0 <- tu %>% distinct() %>%
        mutate(date = as.Date(format(timestamp, '%Y-%m-%d'))) %>%
        group_by(date, trip_id, vehicle_id) %>%
        do((.) %>% trip_times_to_travel_times())
    save(tu0, file = "tu0.rda")
}

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

tu1 <- tu0 %>% mutate(
    short_route_id = gsub('-.*', '', route_id),
    route_number = rnums[route_id]
)

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

rids <- names(sort(table(tu1$short_route_id), TRUE))
RID <- rids[21]; routei %>% filter(grepl(RID, route_id))

tu1r <- tu1 %>% 
    filter(
        short_route_id == RID &
        between(travel_time, 20, 60*30)
    ) %>%
    filter(
        segment_index < max((.)$segment_index)
    ) %>%
    group_by(segment_index) %>%
    do((.) %>% mutate(
        time = as.POSIXct(paste(Sys.Date(), format(time_from, '%H:%M:%S'))),
        dow = as.numeric(format(time_from, '%u')),
        week = as.factor(ifelse(as.numeric(format(time_from, '%d')) < 8, 1, 2)),
        weekend = as.numeric(dow > 5),
        travel_time.std = std(travel_time)
    )) %>%
    ungroup() %>%
    arrange(segment_index, time_from)
 
ggplot(
    tu1r, 
    aes(time, travel_time/60, colour = as.factor(dow))
) + geom_point() + facet_grid(segment_index~weekend, scales = "free_y") +
    scale_x_datetime(label = function(x) format(x, '%H')) +
    geom_smooth(se = FALSE)

ggplot(
    tu1r %>% filter(segment_index == 33), 
    aes(time, travel_time/60, colour = as.factor(dow))
) + geom_point() + facet_grid(week~weekend) +
    scale_x_datetime(label = function(x) format(x, '%H'))

ggplot(
    tu1r %>% filter(segment_index == 33), 
    aes(time, travel_time.std/60, colour = as.factor(dow))
) + geom_point() + facet_grid(week~weekend) +
    scale_x_datetime(label = function(x) format(x, '%H'))

## num -> [30min intervals]
times <- pretty(tu1r$time, 20 * 2)
smry <- tu1r %>%
    mutate(
        period = cut((.)$time, breaks = times)
    ) %>%
    group_by(date, segment_index, period) %>%
    summarize(
        tt = mean(travel_time, na.rm = TRUE),
        sd.tt = sd(travel_time, na.rm = TRUE),
        n = n(),
        dow = first(dow),
        week = first(week),
        weekend = first(weekend)
    ) %>%
    mutate(
        time = as.POSIXct(period)
    )

ggplot(smry %>% filter(segment_index == 18), 
    aes(time, tt/60, colour = as.factor(dow))) + 
    geom_point() +
    facet_grid(week~weekend) +
    scale_x_datetime(label = function(x) format(x, '%H')) 


## drop it into a matrix




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

## short trip ...
tu1 <- tu1 %>% mutate(short_trip_id = gsub('-.*', '', trip_id))
tsegs <- tsegs %>% mutate(short_trip_id = gsub('-.*', '', trip_id))
sx <- tu1 %>% inner_join(tsegs, by = c("short_trip_id", "segment_index" = "shape_road_sequence")) %>%
    filter(travel_time > 20 & travel_time < 60*30) %>%
    arrange(road_segment_id, time_from) %>%
    mutate(speed = length / travel_time) %>%
    filter(speed < 30)

NX1segs <- sx %>% ungroup() %>% 
    filter(grepl(RID, route_id)) %>%
    select(road_segment_id, segment_index) %>% distinct() %>%
    arrange(segment_index) %>%
    filter(segment_index < max(segment_index))

# seg20 <- sx$road_segment_id %>% table %>% sort(TRUE) %>% names
sxnex <- sx %>% filter(road_segment_id %in% NX1segs$road_segment_id) %>%
    mutate(
        segment = factor(road_segment_id, levels = NX1segs$road_segment_id),
        time = as.POSIXct(paste(Sys.Date(), format(time_from, '%H:%M:%S'))),
        dow = as.numeric(format(time_from, '%u')),
        week = as.factor(ifelse(as.numeric(format(time_from, '%d')) < 8, 1, 2)),
        weekend = as.numeric(dow > 5)
    )

times2 <- pretty(sxnex$time, 20 * 2)
smry2 <- sxnex %>% ungroup () %>%
    mutate(
        period = cut((.)$time, breaks = times2)
    ) %>%
    group_by(date, segment, period) %>%
    summarize(
        tt = mean(travel_time, na.rm = TRUE),
        sd.tt = sd(travel_time, na.rm = TRUE),
        n = n(),
        speed = mean(speed, na.rm = TRUE)
    ) %>%
    ungroup () %>%
    mutate(
        time = as.POSIXct(period),
        dow = as.numeric(format(date, '%u')),
        week = as.factor(ifelse(as.numeric(format(date, '%d')) < 8, 1, 2)),
        weekend = as.numeric(dow > 5)
    )


ggplot(smry2, aes(time, speed / 1000 * 60 * 60)) + 
# ggplot(smry2, aes(time, tt / 60)) + 
    geom_point(aes(colour = as.factor(week))) +
    # geom_path(aes(colour = date)) +
    facet_grid(segment~weekend, scales = "free_y") +
    xlab("Time") + 
    ylab("Travel time (min)") +
    scale_x_datetime(label = function(x) format(x, '%H:%M'))
    # geom_smooth(aes(colour = route_number), se = FALSE) + 
    # scale_colour_discrete(guide = FALSE) + 
    # scale_y_log10()
    # geom_smooth(method = lm, formula = y ~ splines::bs(x, 4))


### I want the model to work on both mu and sigma
# i.e., mu(t) and sigma(t)

alldata <- smry2 %>% 
    select(date, time, segment, tt, speed, dow, week, weekend) %>%
    mutate(
        date = factor(date),
        dow = factor(dow),
        dayofweek = factor(format(as.Date(date), '%a'), 
            levels = c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')),
        week = factor(week),
        saturday = as.numeric(dow == 6),
        sunday = as.numeric(dow == 7)
    ) #%>% filter(segment == '4691')
train <- alldata %>% filter(week == 1)
test <- alldata %>% filter(week == 2)

rmse <- function(fit, data) 
    mean((predict(fit, newdata = data) - data$tt)^2, na.rm = TRUE)


fit0 <- lm(tt ~ 1, data = train)
rmse(fit0, test)

fit1 <- lm(tt ~ segment * dow - 1, data = train)
rmse(fit1, test)

fit2 <- lm(tt ~ segment + segment*saturday + segment*sunday - 1, data = train)
rmse(fit2, test)


## spline models
library(splines)
fitSpline <- function(i, j)
    lm(speed ~ bs(time, intercept = TRUE, degree = i):segment + 
        bs(time, intercept = TRUE, degree = j):weekend:segment - 1, data = train)

RMSE <- expand.grid(i = 1:20, j = 1:15)
RMSE$rmse <- apply(RMSE, 1, function(x) {
    rmse(fitSpline(x["i"], x["j"]), test)
})
ggplot(RMSE, aes(i, j, size = log(rmse))) + geom_point()

RMSE[which.min(RMSE$rmse), ]
sfit <- fitSpline(5, 2)
summary(sfit)
test %>% mutate(speed.pred = predict(sfit, newdata = test)) %>%
    ggplot(aes(time, speed.pred*3.6)) + 
    geom_point(aes(y = speed*3.6), data = train, col = 'gray', size = 0.5) +
    geom_point(aes(y = speed*3.6), data = test, col = 'orangered', size = 0.5) +
    geom_path(aes(group = date), colour = 'steelblue') +
    facet_grid(segment~dayofweek, scales = "free_y") +
    xlab("Time") + 
    ylab("Speed (m/s)") + #ylim(0, 100) +
    scale_x_datetime(label = function(x) format(x, '%H:%M'))
    # ylab("Travel time (min)") +


## lag model
alldata <- alldata %>% mutate(time_index = time %>% as.factor %>% as.integer)

alldata.lag <- alldata %>%
    filter(time_index > 1) %>%
    mutate(time_index_lag = time_index - 1) %>%
    left_join(
        alldata %>% select(date, segment, time_index, tt, speed),
        by = c('date', 'segment', 'time_index_lag' = 'time_index'),
        suffix = c("", ".lag")
    )


train.lag <- alldata.lag %>% filter(week == 1)
test.lag <- alldata.lag %>% filter(week == 2)

ggplot(train.lag, aes(time, tt - tt.lag)) +
    geom_point() +
    facet_grid(segment~dayofweek, scales = "free_y") +
    scale_x_datetime(label = function(x) format(x, '%H:%M'))


fitLagSpline <- function(i, j)
    lm(I(tt - tt.lag) ~ bs(time, degree = i, intercept = TRUE):segment +
        bs(time, degree = j, intercept = TRUE):weekend:segment - 1, data = train.lag)

RMSE <- expand.grid(i = 1:10, j = 1:5)
RMSE$rmse <- apply(RMSE, 1, function(x) {
    rmse(fitLagSpline(x["i"], x["j"]), test.lag)
})
ggplot(RMSE, aes(i, j, size = rmse)) + geom_point()

RMSE[which.min(RMSE$rmse), ]
lfit <- fitLagSpline(8, 3)
summary(lfit)

test.lag %>% mutate(tt.pred = predict(lfit, newdata = test.lag)) %>%
    ggplot(aes(time, (tt.pred)/60)) + 
    geom_point(aes(y = (tt - tt.lag)/60), data = train.lag, col = 'gray', size = 0.5) +
    geom_path(aes(y = (tt.pred)/60, group = date), 
        data = train.lag %>% mutate(tt.pred = predict(lfit, train.lag)), col = 'gray') +
    geom_point(aes(y = (tt - tt.lag)/60), col = 'orangered', size = 0.5) +
    # geom_path(col = 'orangered') +
    facet_grid(segment~weekend, scales = "free_y") +
    xlab("Time") + 
    ylab("Travel Time (min)") + 
    scale_x_datetime(label = function(x) format(x, '%H:%M'))
    # ggplot(aes(time, speed.pred*3.6)) + 
    # geom_point(aes(y = speed*3.6), data = train, col = 'gray', size = 0.5) +
    # geom_point(aes(y = speed*3.6), data = test, col = 'orangered', size = 0.5) +
    # ylab("Travel time (min)") +


predict(fit, test.lag %>% mutate(tt.lag = 0))