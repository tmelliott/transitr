## Read ETAs from pb files in simulations/simXXX/etas
## using a cpp function to speed things up ... 
# library(devtools)
library(tidyverse)

ts2dt <- function(ts) as.POSIXct(ts, origin = "1970-01-01")

# Load the actual arrival times?
tudir <- file.path("simulations", "archive")
tufiles <- list.files(tudir, pattern = "^trip_updates.+\\.pb$", full.names = TRUE)
# system.time( transitr:::processEtas(tufiles, "arrivaldata.csv", "at_gtfs.db") )

system.time(
    arrivaldata <- readr::read_csv(
        "arrivaldata.csv",
        col_names = c(
            "trip_id", "route_id", "vehicle_id", "timestamp", 
            "stop_sequence", "current_delay", 
            "arrival_time", "scheduled_arrival",
            "lower", "upper", "type"
        ),
        col_types = "ccciiiiiiif"
    ) %>% 
    mutate(
        scheduled_arrival = ts2dt(scheduled_arrival),
        delay = current_delay,
        arrival_time = scheduled_arrival + delay
    ) %>%
    select(trip_id, route_id, vehicle_id, stop_sequence, 
        scheduled_arrival, arrival_time, delay) %>%
    unique() %>%
    group_by(trip_id) %>%
    do({
        x <- (.) %>% unique() %>% 
            filter(between(delay, -30*60, 60*60))
        x <- x %>% group_by(stop_sequence) %>%
            do({
                ## the smallest absolute delay
                z <- (.)
                z[which.min(abs(z$delay)), ]
            })
        vids <- table(x$vehicle_id)
        x %>% filter(vehicle_id == names(vids)[which.max(vids)])
    })
)


tids <- arrivaldata %>% group_by(trip_id) %>%
    summarize(n_vehicle = length(unique(vehicle_id))) %>%
    filter(n_vehicle > 1) %>%
    pull(trip_id)

## quick inspection of the raw arrival data
# egg::ggarrange(
# ggplot(arrivaldata %>% filter(!trip_id %in% tids) %>% arrange(arrival_time), 
#     aes(arrival_time, vehicle_id, group = trip_id, colour = trip_id)) +
#     xlab("time") + theme(legend.position = "none") +
#     geom_path(),
# ggplot(arrivaldata %>% filter(trip_id %in% tids) %>% arrange(arrival_time), 
#     aes(arrival_time, vehicle_id, group = trip_id, colour = trip_id)) +
#     xlab("time") + theme(legend.position = "none") +
#     geom_path(),
# ncol = 2
# )


sim <- "sim000"
etadir <- file.path("simulations", sim, "etas")
etafiles <- list.files(etadir, pattern = ".pb$", full.names = TRUE)

# load_all()
system.time( transitr:::processEtas(etafiles, "etas.csv", "at_gtfs.db") )

system.time(
    etas <- readr::read_csv(
        "etas.csv",
        col_names = c(
            "trip_id", "route_id", "vehicle_id", "timestamp", 
            "stop_sequence", "current_delay", 
            "arrival_time", "scheduled_arrival",
            "lower", "upper", "type"
        ),
        col_types = "ccciiiiiiif"
    )
)
etas <- etas %>% 
    mutate(
        timestamp = ts2dt(timestamp),
        eta_prediction = ts2dt(arrival_time),
        scheduled_arrival = ts2dt(scheduled_arrival),
        lower_prediction = ts2dt(lower),
        upper_prediction = ts2dt(upper)
        # gtfs_eta = scheduled_arrival + current_delay
    ) %>%
    select(trip_id, route_id, timestamp, stop_sequence, eta_prediction,
        current_delay, lower_prediction, upper_prediction)


## silly plot of the ETAs
# for (TRIP_ID in unique(etas$trip_id)) {
#     ad <- arrivaldata %>% filter(trip_id == TRIP_ID)
#     p <- ggplot(etas %>% filter(trip_id == TRIP_ID)) +
#         geom_point(aes(eta_prediction, timestamp)) +
#         geom_segment(
#             aes(
#                 x = lower_prediction, y = timestamp,
#                 xend = upper_prediction, yend = timestamp
#             )
#         ) +
#         geom_vline(aes(xintercept = arrival_time),
#             data = arrivaldata %>% filter(trip_id == TRIP_ID),
#             colour = "orangered"
#         ) +
#         geom_vline(aes(xintercept = scheduled_arrival),
#             data = ad
#         ) +
#         ggtitle(TRIP_ID) +
#         facet_wrap(~stop_sequence)
#     print(p)
#     grid::grid.locator()
# }


library(RSQLite())
con <- dbConnect(SQLite(), "at_gtfs.db")
HMS2sec <- function(hms) {
    sapply(strsplit(hms, ":"), function(z) {
        sum(as.integer(z) * c(60*60, 60, 1))
    })
}
t00 <- paste(format(arrivaldata$arrival_time[1], "%Y-%m-%d"), "00:00:00") %>%
    as.POSIXct %>% as.integer
adata <- arrivaldata %>% select(trip_id, stop_sequence, scheduled_arrival, arrival_time, delay)

# ## convert this to a form [trip_id, timestamp, delay]
# adata %>% full_join(etas, by = c("trip_id"))

#  group_by(trip_id, arrival_time) %>%
#     group_modify(~{
#         tid <- .y$trip_id
#         ts <- etas %>% filter(trip_id == tid) %>% pull(timestamp) %>% unique()
#         if (length(ts) == 0) return(NULL)
#         .x
#     })

eta_data <- etas %>% select(-current_delay) %>% 
    # filter(trip_id %in% unique(etas$trip_id)) %>%
    left_join(adata %>% rename(actual_arrival = arrival_time) %>%
            select(trip_id, stop_sequence, scheduled_arrival, actual_arrival, delay),
        by = c("trip_id", "stop_sequence")
    ) %>%
    mutate(
        time_until_arrival = as.integer(actual_arrival - timestamp),
        eta = as.integer(eta_prediction - timestamp)
    ) %>%
    group_by(trip_id) %>%
    group_modify(~ {
        # .y <- tibble(trip_id = unique(etas$trip_id)[2])
        # .x <- eta_data %>% filter(trip_id == .y$trip_id)

        d <- adata %>% filter(trip_id == .y$trip_id & !is.na(delay))
        dat <- tibble(timestamp = unique(.x$timestamp)) %>%
            mutate(
                current_delay = sapply(timestamp, 
                    function(t) {
                        suppressWarnings(a <- d$delay[max(which(d$arrival_time <= t))])
                        ifelse(is.na(a), 0, a)
                    }
                )
            )

        .x %>% left_join(dat, by = "timestamp") %>% arrange(stop_sequence) %>% 
            mutate(
                gtfs_arrival_time = scheduled_arrival + current_delay,
                gtfs_eta = as.integer(scheduled_arrival + current_delay) - as.integer(timestamp)
            )
    }) %>% 
    filter(time_until_arrival > 0)

## A VERY VERY BASIC RMSE THING
RMSE <- list(
    transitr = sqrt(mean(eta_data$eta^2, na.rm = TRUE)),
    ci_cov = mean(with(eta_data, actual_arrival > lower_prediction & actual_arrival < upper_prediction), na.rm = TRUE),
    gtfs = sqrt(mean(eta_data$gtfs_eta^2, na.rm = TRUE))
)
RMSE

range(eta_data$timestamp)

# q("no")
for (TRIP in unique(eta_data$trip_id)) {
    print(TRIP)
    routedata <- eta_data %>% 
        filter(trip_id == TRIP & !is.na(eta))# & timestamp > as.POSIXct("2019-08-19 10:00:00"))
    if (nrow(routedata) == 0) next()
    p <- ggplot(routedata %>% filter(time_until_arrival > 0), 
        aes(timestamp)) +
        geom_pointrange(aes(y = eta_prediction, ymin = lower_prediction, ymax = upper_prediction)) +
        geom_point(aes(y = gtfs_arrival_time), colour = "red") +
        geom_hline(aes(yintercept = actual_arrival), colour = "blue") +
        geom_hline(aes(yintercept = scheduled_arrival), colour = "black", lty = 2) +
        facet_wrap(~stop_sequence,scales="free_y") +
        theme(legend.position = "none")
    print(p)
    grid::grid.locator()
}


## RAW eta state data
raw <- readr::read_csv(
        "simulations/sim000/eta_state.csv",
        col_names = c(
            "trip_id", "stop_sequence", "timestamp",
            "X", "P", "Xhat", "Phat", "Z", "E", "Xnew", "Pnew"
        ),
        col_types = "ciinnnnnnnn"
    ) %>% filter(X < 2*60*60) %>%
    mutate(
        timestamp = ts2dt(timestamp),
        stop_sequence = stop_sequence + 1
    )
max(raw$timestamp)


get_ci <- function(x, p, which = c("lower", "upper"), q = 0.95, ts) {
    which <- match.arg(which)
    q <- (1 - q) / 2
    if (which == "upper") q <- 1 - q
    p <- ifelse(p > 1e6, 0, p)
    z <- qnorm(q, x, sqrt(p))
    if (missing(ts)) return(z)
    t <- as.POSIXct(ts + z, origin = "1970-01-01")
    if (which == "upper") t <- t + 60
    as.POSIXct(
        paste(
            format(t, "%Y-%m-%d %H:"),
            as.integer(format(t, "%M")),
            ":00"
        )
    )
}

## pick a trip
# tid <- "1141155994-20190806160740_v82.21"
# tid <- "1141156213-20190806160740_v82.21"
# tid <- "1141156521-20190806160740_v82.21"

for (tid in unique(raw$trip_id)) {
    trip_data <- raw %>% filter(trip_id == tid)
    arr_data <- arrivaldata %>% filter(trip_id == tid)

    p <- ggplot(trip_data, aes(timestamp)) +
        geom_hline(aes(x = NULL, y = NULL, yintercept = arrival_time),
            data = arr_data) +
        geom_vline(aes(x = NULL, y = NULL, xintercept = arrival_time),
            data = arr_data) +
        geom_hline(aes(x = NULL, y = NULL, yintercept = scheduled_arrival),
            data = arr_data, lty = 3) +
        geom_linerange(
            aes(
                x = timestamp,
                # y = timestamp + Z,
                ymin = timestamp + get_ci(Z, E, "lower"),
                ymax = timestamp + get_ci(Z, E, "upper")
            ), 
            colour = "red"
        ) +
        geom_linerange(
            aes(
                x = timestamp+5,
                ymin = timestamp + get_ci(Xhat, Phat, "lower"),
                ymax = timestamp + get_ci(Xhat, Phat, "upper")
            ), 
            colour = "blue"
        ) +
        geom_linerange(
            aes(
                x = timestamp+10,
                # y = NULL,#timestamp + Xnew,
                ymin = get_ci(Xnew, Pnew, "lower", ts = timestamp),
                ymax = get_ci(Xnew, Pnew, "upper", ts = timestamp)
            )
        ) +
        facet_wrap(~stop_sequence, scale = "free_y")
    plot(p)
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


