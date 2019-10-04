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
            "arrival_time", "scheduled_arrival"
        ),
        col_types = "cciiii"
    )
)
ts2dt <- function(ts) as.POSIXct(ts, origin = "1970-01-01")
etas <- etas %>% 
    mutate(
        timestamp = ts2dt(timestamp),
        arrival_time = ts2dt(arrival_time),
        scheduled_arrival = ts2dt(scheduled_arrival)
    )

load("simulations/arrivaldata.rda")
arrivaldata <- arrivaldata %>%
    ungroup() %>%
    spread(key = type, value = time) %>%
    mutate(
        actual_arrival = ts2dt(ifelse(is.na(arrival), departure, arrival))
    )

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
    filter(between(delay, -60*60, 60*60)) %>%
    filter(eta < 2*60*60 & time_until_arrival >= 0)

ggplot(eta_data, aes(time_until_arrival, (eta - actual_arrival))) +
    geom_point() +
    # facet_wrap(~route_id, scales = "free") +
    geom_abline(colour = "orangered")

bad_trip <- eta_data %>% 
    filter(trip_id=="51446171528-20190806160740_v82.21") %>% 
    arrange(timestamp, stop_sequence)

ggplot(bad_trip, aes(stop_sequence, arrival_time)) +
    geom_point()+
    facet_wrap(~timestamp)


