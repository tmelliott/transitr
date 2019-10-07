## Read ETAs from pb files in simulations/simXXX/etas
## using a cpp function to speed things up ... 
# library(devtools)
library(tidyverse)

sim <- "sim000"
etadir <- file.path("simulations", sim, "etas")
etafiles <- list.files(etadir, pattern = ".pb$", full.names = TRUE)

# load_all()
system.time( transitr:::processEtas(etafiles, "etas.csv", "at_gtfs.db") )

system.time(
    etas <- readr::read_csv(
        "etas.csv",
        col_names = c(
            "trip_id", "route_id", "timestamp", "stop_sequence", 
            "current_delay", "arrival_time", "scheduled_arrival"
        ),
        col_types = "cciiiii"
    )
)
ts2dt <- function(ts) as.POSIXct(ts, origin = "1970-01-01")
etas <- etas %>% 
    mutate(
        timestamp = ts2dt(timestamp),
        eta_prediction = ts2dt(arrival_time),
        scheduled_arrival = ts2dt(scheduled_arrival)
        # gtfs_eta = scheduled_arrival + current_delay
    ) %>%
    select(trip_id, route_id, timestamp, stop_sequence, eta_prediction,
        current_delay)

## Load the actual arrival times?
# tudir <- file.path("simulations", "archive")
# tufiles <- list.files(tudir, pattern = "^trip_updates.+\\.pb$", full.names = TRUE)
# system.time( transitr:::processEtas(tufiles, "arrivaldata.csv", "at_gtfs.db") )

# system.time(
#     arrivaldata <- readr::read_csv(
#         "arrivaldata.csv",
#         col_names = c(
#             "trip_id", "route_id", "timestamp", "stop_sequence", 
#             "current_delay", "arrival_time", "scheduled_arrival"
#         ),
#         col_types = "cciiiii"
#     )
# )
# arrivaldata <- arrivaldata %>% 
#     mutate(
#         scheduled_arrival = ts2dt(scheduled_arrival),
#         delay = current_delay,
#         arrival_time = scheduled_arrival + delay
#     ) %>%
#     select(trip_id, route_id, stop_sequence, scheduled_arrival, arrival_time, delay) %>%
#     unique()


## silly plot of the ETAs
for (TRIP_ID in unique(etas$trip_id)) {
    p <- ggplot(etas %>% filter(trip_id == TRIP_ID)) +
        geom_point(aes(eta_prediction, timestamp)) +
        geom_vline(aes(xintercept = arrival_time),
            data = arrivaldata %>% filter(trip_id == TRIP_ID),
            colour = "orangered"
        ) +
        geom_vline(aes(xintercept = scheduled_arrival),
            data = arrivaldata %>% filter(trip_id == TRIP_ID)
        ) +
        ggtitle(TRIP_ID) +
        facet_wrap(~stop_sequence)
    print(p)
    grid::grid.locator()
}


library(RSQLite())
con <- dbConnect(SQLite(), "at_gtfs.db")
HMS2sec <- function(hms) {
    sapply(strsplit(hms, ":"), function(z) {
        sum(as.integer(z) * c(60*60, 60, 1))
    })
}
t00 <- paste(format(arrivaldata$arrival_time[1], "%Y-%m-%d"), "00:00:00") %>%
    as.POSIXct %>% as.integer
# adata <- arrivaldata %>%
#     group_by(trip_id) %>%
#     do({
#         x <- (.) %>% rename(actual_arrival = arrival_time)
#         tid <- x$trip_id[1]
#         Schedule <- con %>% tbl("stop_times") %>%
#             filter(trip_id == !!tid) %>%
#             select(stop_sequence, arrival_time) %>%
#             arrange(stop_sequence) %>% collect()
#         x2 <- left_join(Schedule, 
#                 x %>% select(stop_sequence, actual_arrival),
#                 by = "stop_sequence"
#             ) %>%
#             mutate(
#                 arrival_time = ts2dt(HMS2sec(arrival_time) + t00),
#                 arrival_delay = as.integer(actual_arrival) - as.integer(arrival_time)
#             )
#         if (is.na(x2$arrival_delay[1])) x2$arrival_delay[1] <- 0
#         for (i in 2:nrow(x2)) {
#             if (is.na(x2$arrival_delay[i])) 
#                 x2$arrival_delay[i] <- x2$arrival_delay[i-1]
#         }
#         x2
#     })
adata <- arrivaldata %>% select(trip_id, stop_sequence, scheduled_arrival, arrival_time, delay)

eta_data <- etas %>% select(-current_delay) %>%
    left_join(adata %>% rename(actual_arrival = arrival_time) %>%
            select(trip_id, stop_sequence, scheduled_arrival, actual_arrival, delay),
        by = c("trip_id", "stop_sequence")
    ) %>%
    mutate(
        time_until_arrival = as.integer(actual_arrival - timestamp),
        eta = as.integer(eta_prediction - timestamp)
    ) %>%
    group_by(trip_id, timestamp) %>%
    do({
        x <- (.) %>% arrange(stop_sequence)
        # x <- eta_data %>% 
        #     filter(trip_id=="8301168082-20190806160740_v82.21" & timestamp == "2019-08-19 05:00:11") %>% 
        #     arrange(stop_sequence)
        curdelay <- adata %>% 
            filter(trip_id == x$trip_id[1] & arrival_time <= x$timestamp[1]) %>%
            filter(!is.na(delay))
        if (nrow(curdelay) == 0) curdelay <- 0
        else curdelay <- curdelay$delay[which.max(curdelay$arrival_time)]
        x %>% mutate(
            current_delay = curdelay, 
            gtfs_arrival_time = scheduled_arrival + current_delay,
            gtfs_eta = as.integer(scheduled_arrival + current_delay) - as.integer(timestamp)
        )
    })

eta_data <- eta_data %>% filter(time_until_arrival > 0)

    #     gtfs_eta = as.integer(scheduled_arrival + arrival_delay - timestamp)
    # ) %>%
    # filter(between(delay, -60*60, 60*60) & time_until_arrival > 0)# & current_delay != 0) # %>%
    # # filter(eta < 2*60*60 & between(time_until_arrival, 0, 60*60*2))

## A VERY VERY BASIC RMSE THING
RMSE <- list(
    transitr = sqrt(mean(eta_data$eta^2)),
    gtfs = sqrt(mean(eta_data$gtfs_eta^2, na.rm = TRUE))
)


for (TRIP in unique(eta_data$trip_id)) {
    routedata <- eta_data %>% filter(trip_id == TRIP)

    p <- ggplot(routedata %>% filter(time_until_arrival > 0), 
        aes(timestamp)) +
        geom_point(aes(y = eta_prediction)) +
        geom_point(aes(y = gtfs_arrival_time), colour = "red") +
        geom_hline(aes(yintercept = actual_arrival), colour = "blue") +
        geom_hline(aes(yintercept = scheduled_arrival), colour = "black", lty = 2) +
        facet_wrap(~stop_sequence,scales="free_y") +
        theme(legend.position = "none")
    print(p)
    grid::grid.locator()
}

    # p <- ggplot(routedata, aes(time_until_arrival/60, 
    #     as.integer(arrival_time - actual_arrival)/60)) +
    #     geom_point(aes(colour = trip_id)) +
    #     geom_point(aes(y = as.integer(gtfs_eta - actual_arrival)/60),
    #         colour = "red", size = 1) +
    #     facet_wrap(~stop_sequence)
    #     # geom_abline(colour = "orangered")
    # print(p)
    # grid::grid.locator()

    # p <- ggplot(routedata, aes(timestamp, 
    #     as.integer(arrival_time-timestamp)/60)) +
    #     geom_point(aes(colour = as.factor(stop_sequence))) +
    #     geom_point(aes(y = as.integer(gtfs_eta - timestamp)/60), colour = "red", size = 0.5)+
    #     facet_wrap(~stop_sequence)
    # print(p)
    grid::grid.locator()
}

bad_trip <- eta_data %>% 
    filter(trip_id=="51357146239-20190806160740_v82.21") %>% 
    # filter(trip_id=="51446171528-20190806160740_v82.21") %>% 
    arrange(timestamp, stop_sequence)

## segment ids for this route
## funs from view_network.R
segs <- get_segment_data(route_ids = bad_trip$route_id[1])
dat <- get_segments("simulations/sim000/segment_states.csv") %>%
    filter(segment_id %in% segs$road_segment_id) %>%
    left_join(segs, by = c("segment_id" = "road_segment_id"))

ggplot(dat, aes(timestamp, travel_time)) +
    geom_pointrange(aes(
        ymin = travel_time - sqrt(uncertainty),
        ymax = travel_time + sqrt(uncertainty)
    )) +
    facet_wrap(~segment_id, scales = "free_y")

ggplot(bad_trip, aes(timestamp, arrival_time)) +
    geom_point()+
    facet_wrap(~stop_sequence)


