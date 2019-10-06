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
        arrival_time = ts2dt(arrival_time),
        scheduled_arrival = ts2dt(scheduled_arrival),
        gtfs_eta = scheduled_arrival + current_delay
    )

load("simulations/arrivaldata.rda")
arrivaldata <- arrivaldata %>%
    ungroup() %>%
    spread(key = type, value = time) %>%
    mutate(
        actual_arrival = ts2dt(ifelse(is.na(arrival), departure, arrival))
    )# %>% 
    ## Add AT estimate of arrival time 
# arrivaldata %>%
#     group_by(trip_id) %>%
#     do({
#         x <- (.)
#         tid <- x$trip_id[1]
#         Schedule <- con %>% tbl("stop_times") %>%
#             filter(trip_id == !!tid) %>%
#             select(stop_sequence, arrival_time) %>%
#             arrange(stop_sequence) %>% collect()
#         x2 <- left_join(Schedule, x %>% select(stop_sequence, actual_arrival)) %>%
#             mutate(
#                 arrival_time = 
#                     as.POSIXct(paste(
#                         format(actual_arrival[1], "%Y-%m-%d"), 
#                         arrival_time
#                     )),
#                 arrival_delay = as.integer(actual_arrival - arrival_time)
#             )
#         x2
#     })

eta_data <- etas %>% 
    left_join(
        arrivaldata %>% select(trip_id, stop_sequence, actual_arrival),
        by = c("trip_id", "stop_sequence")
    ) %>%
    mutate(
        time_until_arrival = as.integer(actual_arrival - timestamp),
        eta = as.integer(arrival_time - timestamp),
        delay = as.integer(actual_arrival - scheduled_arrival)
    ) %>%
    filter(between(delay, -60*60, 60*60) & time_until_arrival > 0)# & current_delay != 0) # %>%
    # filter(eta < 2*60*60 & between(time_until_arrival, 0, 60*60*2))

## A VERY VERY BASIC RMSE THING
RMSE_ME <- sqrt(mean(with(eta_data, as.integer(arrival_time - actual_arrival))^2))
RMSE_GTFS <- sqrt(mean(with(eta_data, as.integer(gtfs_eta - actual_arrival))^2))



ggplot(eta_data, aes(time_until_arrival/60, 
    as.integer(arrival_time - actual_arrival)/60)) +
    geom_point(aes(colour = stop_sequence)) +
    geom_point(aes(y = as.integer(gtfs_eta - actual_arrival)/60),
        colour = "red", size = 1) +
    facet_wrap(~route_id, scales = "free") 
    # geom_abline(colour = "orangered")

ggplot(eta_data, aes(timestamp, 
    as.integer(arrival_time-timestamp)/60)) +
    geom_point(aes(colour = as.factor(stop_sequence))) +
    geom_point(aes(y = as.integer(gtfs_eta - timestamp)/60), colour = "red", size = 0.5)+
    facet_wrap(~route_id)


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


