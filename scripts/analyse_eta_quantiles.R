library(tidyverse)
ts2dt <- function(ts) as.POSIXct(ts, origin = "1970-01-01")
load("simulations/raw_arrivaldata.rda")
rawarrival <- arrivaldata
load("simulations/arrivaldata_etas.rda")

eta_quantiles <-
    read_csv("simulations/sim002/eta_quantiles.csv",
        col_names = c("trip_id", "vehicle_id", "stop_sequence", "timestamp", "eta", "quantile"),
        col_types = "cciiin"
    ) %>%
    left_join(
        arrivaldata %>%
            rename(actual_arrival = arrival_time) %>%
            select(trip_id, vehicle_id, stop_sequence, scheduled_arrival, actual_arrival, delay),
        by = c("trip_id", "vehicle_id", "stop_sequence")
    ) %>%
    mutate(
        scheduled_arrival = as.integer(scheduled_arrival),
        actual_arrival = as.integer(actual_arrival),
        time_until_arrival = actual_arrival - timestamp
    ) %>%
    group_by(trip_id, vehicle_id) %>%
    group_modify(~ {
        # .y <- tibble(trip_id = unique(etas$trip_id)[2])
        # .x <- eta_data %>% filter(trip_id == .y$trip_id)

        d <- rawarrival %>%
            filter(
                trip_id == .y$trip_id &
                vehicle_id == .y$vehicle_id &
                !is.na(delay)
            )
        dat <- tibble(timestamp = unique(.x$timestamp)) %>%
            mutate(
                current_delay = sapply(
                    timestamp,
                    function(t) {
                        suppressWarnings(a <- d$delay[max(which(d$arrival_time <= t))])
                        ifelse(is.na(a), 0, a)
                    }
                )
            )

        .x %>%
            left_join(dat, by = "timestamp") %>%
            mutate(
                scheduled_arrival = as.integer(scheduled_arrival),
                timestamp = as.integer(timestamp)
            ) %>%
            arrange(stop_sequence) %>%
            mutate(
                gtfs_arrival_time = scheduled_arrival + current_delay,
                gtfs_eta = scheduled_arrival + current_delay - timestamp
            )
    }) %>%
    filter(time_until_arrival >= 0) %>%
    mutate(timestamp = ts2dt(timestamp))

eta_smry <- eta_quantiles %>% ungroup %>%
    group_by(trip_id, vehicle_id, stop_sequence, timestamp) %>%
    group_modify(~{
        etaq50 <- with(.x, eta[max(which(quantile <= 0.5))])
        etaq40 <- with(.x, eta[max(which(quantile <= 0.4))])
        etaq25 <- with(.x, eta[max(which(quantile <= 0.25))])
        tibble(
            err50 = .x$time_until_arrival[1] / 60 - etaq50,
            err40 = .x$time_until_arrival[1] / 60 - etaq40,
            err25 = .x$time_until_arrival[1] / 60 - etaq25,
            err_gtfs = (.x$time_until_arrival[1] - .x$gtfs_eta[1]) / 60
        )
    })

ggplot(eta_smry, aes(stop_sequence, err50)) + geom_point()

i <- 10000
T <- eta_quantiles$trip_id[i]
S <- eta_quantiles$stop_sequence[i]
TS <- eta_quantiles$timestamp[i]

eta1 <- eta_quantiles %>%
    filter(trip_id == T & stop_sequence == S & timestamp == TS) %>%
    mutate(p = c(0, diff(quantile)))

etaq50 <- with(eta1, eta[max(which(quantile <= 0.5))])
etaq40 <- with(eta1, eta[max(which(quantile <= 0.4))])
etaq25 <- with(eta1, eta[max(which(quantile <= 0.25))])

ggplot(eta1, aes(eta, quantile)) +
    geom_path() +
    geom_point() +
    geom_hline(yintercept = c(0.25, 0.4, 0.5), lty = 2) +
    geom_vline(xintercept = c(etaq25, etaq40, etaq50), lty = 3) +
    geom_vline(xintercept = eta1$time_until_arrival[1] / 60, colour = "orangered") +
    geom_vline(xintercept = eta1$gtfs_eta[1] / 60, colour = "magenta")
