### Load using 'data.table'
library(data.table)

eta_quantiles <- fread(
    "simulations/sim_100/eta_quantiles.csv",
    colClasses = c(
        "factor", "factor", "integer", "integer",
        "integer", "numeric"),
    col.names = c(
        "trip_id", "vehicle_id", "stop_sequence", "timestamp",
        "eta", "quantile"
    )
)

## some more setup
library(tidyverse)
ts2dt <- function(ts) as.POSIXct(ts, origin = "1970-01-01")
load("simulations/raw_arrivaldata.rda")
rawarrival <- arrivaldata %>%
    as_tibble %>%
    mutate(
        trip_id = as.factor(trip_id),
        route_id = as.factor(route_id),
        vehicle_id = as.factor(vehicle_id)
    ) %>%
    as.data.table()
load("simulations/arrivaldata_etas.rda")
arrivaldata <- arrivaldata %>%
    ungroup() %>%
    mutate(
        trip_id = as.factor(trip_id),
        route_id = as.factor(route_id),
        vehicle_id = as.factor(vehicle_id),
        scheduled_arrival = as.integer(scheduled_arrival),
        arrival_time = as.integer(arrival_time)
    ) %>%
    as.data.table()

# gtfs_delays <- arrivaldata %>%
#     group_by(trip_id, vehicle_id) %>%
#     do({
#         # .y <- tibble(trip_id = "1015073057-20190806160740_v82.21", vehicle_id = "7D68")
#         # .x <- arrivaldata %>% filter(trip_id == .y$trip_id & vehicle_id == .y$vehicle_id)

#         .x <- (.)
#         tripid <- as.character(.x$trip_id[1])
#         vehicleid <- as.character(.x$vehicle_id[1])
#         d <- rawarrival %>%
#             filter(
#                 trip_id == tripid &
#                 vehicle_id == vehicleid &
#                 !is.na(delay)
#             )
#         eq <- eta_quantiles[
#             trip_id == tripid & vehicle_id == vehicleid
#         ]
#         tibble(timestamp = unique(eq$timestamp)) %>%
#             mutate(
#                 current_delay = sapply(
#                     timestamp,
#                     function(t) {
#                         suppressWarnings(a <- d$delay[max(which(d$arrival_time <= t))])
#                         ifelse(is.na(a), 0, a)
#                     }
#                 )
#             )
#     })


eta_q_05 <- eta_quantiles[
    quantile < 0.05,
    .(eta_q05 = max(eta)),
    by = .(trip_id, vehicle_id, stop_sequence, timestamp)
]
eta_q_25 <- eta_quantiles[
    quantile < 0.25,
    .(eta_q25 = max(eta)),
    by = .(trip_id, vehicle_id, stop_sequence, timestamp)
]
eta_q_40 <- eta_quantiles[
    quantile < 0.4,
    .(eta_q40 = max(eta)),
    by = .(trip_id, vehicle_id, stop_sequence, timestamp)
]
eta_q_50 <- eta_quantiles[
    quantile < 0.5,
    .(eta_q50 = max(eta)),
    by = .(trip_id, vehicle_id, stop_sequence, timestamp)
]
eta_q_60 <- eta_quantiles[
    quantile < 0.6,
    .(eta_q60 = max(eta)),
    by = .(trip_id, vehicle_id, stop_sequence, timestamp)
]
eta_q_75 <- eta_quantiles[
    quantile < 0.75,
    .(eta_q75 = max(eta)),
    by = .(trip_id, vehicle_id, stop_sequence, timestamp)
]
eta_q_90 <- eta_quantiles[
    quantile < 0.9,
    .(eta_q90 = max(eta)),
    by = .(trip_id, vehicle_id, stop_sequence, timestamp)
]

eta_q <- eta_q_05[
    eta_q_25, on = .(trip_id, vehicle_id, stop_sequence, timestamp)
][
    eta_q_40, on = .(trip_id, vehicle_id, stop_sequence, timestamp)
][
    eta_q_50, on = .(trip_id, vehicle_id, stop_sequence, timestamp)
][
    eta_q_60, on = .(trip_id, vehicle_id, stop_sequence, timestamp)
][
    eta_q_75, on = .(trip_id, vehicle_id, stop_sequence, timestamp)
][
    eta_q_90, on = .(trip_id, vehicle_id, stop_sequence, timestamp)
]

eta_final <- arrivaldata[
    eta_q,
    on = .(trip_id, vehicle_id, stop_sequence)
][!is.na(arrival_time)]

gtfs_times <- arrivaldata[,.(trip_id,vehicle_id,arrival_time,delay)] %>%
    setnames(c("arrival_time", "delay"), c("timestamp", "current_delay"))

gtfs_delay <- integer(nrow(eta_final))
eta_final$gtfs_eta <- 0L
Ntrip <- length(unique(eta_final$trip_id))
ti <- 0
for (trip in unique(eta_final$trip_id)) {
    etrip <- eta_final[trip_id == trip]

    ti <- ti + 1
    Nvehicle <- length(unique(etrip$vehicle_id))
    vi <- 0
    for (vehicle in unique(etrip$vehicle_id)) {
        ev <- etrip[vehicle_id == vehicle]
        gt <- gtfs_times[
            trip_id == trip &
            vehicle_id == vehicle
        ]

        vi <- vi + 1
        Nt <- length(unique(ev$timestamp))
        tti <- 0
        for (t in unique(ev$timestamp)) {
            cur_delay <- gt[timestamp <= t]$current_delay
            eta_final[
                trip_id == trip_id &
                vehicle_id == vehicle &
                timestamp == t,
                gtfs_eta := (scheduled_arrival + cur_delay - timestamp)
            ]

            tti <- tti + 1
            cat(
                glue::glue(
                    "\r Trip {ti}/{Ntrip}; Vehicle {vi}/{Nvehicle}; Time {tti}/{Nt}}"
                )
            )
        }
    }
}

eta_final[
    gtfs_times,
    on = .(trip_id, vehicle_id, timestamp)
][!is.na(arrival_time)]

eta_final[order(trip_id, timestamp, stop_sequence)]


all_trip_ids <- levels(eta_quantiles$trip_id)
for (tid in all_trip_ids) {
    vid <- rawarrival %>%
        filter(trip_id == !!tid) %>%
        group_by(vehicle_id) %>%
        summarize(n = n()) %>%
        filter(n == max(n)) %>%
        pull(vehicle_id) %>%
        as.character()

}


### Load using 'disk.frame'
library(disk.frame)
setup_disk.frame()

fs_dir <- "simulations/sim_100/eta_quantiles.df"
if (dir.exists(fs_dir)) {
    eta_quantiles_df <- disk.frame(fs_dir)
} else {
    eta_quantiles_df <- csv_to_disk.frame(
        "simulations/sim_100/eta_quantiles.csv",
        outdir = fs_dir,
        colClasses = c(
            "character", "character", "integer", "integer",
            "integer", "numeric"),
        col.names = c(
            "trip_id", "vehicle_id", "stop_sequence", "timestamp",
            "eta", "quantile"
        )
    )
}

## some more setup
library(tidyverse)
ts2dt <- function(ts) as.POSIXct(ts, origin = "1970-01-01")
load("simulations/raw_arrivaldata.rda")
rawarrival <- arrivaldata

all_trip_ids <- unique(
    collect(srckeep(eta_quantiles_df, "trip_id"))$trip_id
)

for (trip_id in all_trip_ids) {
    adf <- rawarrival %>%
        filter(trip_id == !!trip_id)

    # vehicle ID:
    vids <- adf %>% group_by(vehicle_id) %>%
        summarize(n = n())
    vid <- vids$vehicle_id[which.max(vids$n)]
    adf <- adf %>% filter(vehicle_id == !!vid)

    # get
    tdf <- eta_quantiles_df %>%
        filter(trip_id == !!trip_id & vehicle_id == !!vid) %>%
        collect()
}


## Now join arrivals data to eta_quantiles:
eta_joined <- eta_quantiles_df %>%
    left_join(
        arrivaldata %>%
            #filter(trip_id == trip) %>%
            rename(actual_arrival = arrival_time) %>%
            select(trip_id, vehicle_id, stop_sequence, scheduled_arrival, actual_arrival, delay),
        by = c("trip_id", "vehicle_id", "stop_sequence"),
        outdir = "simulations/sim_100/eta_quantiles_joined.df"
    )


load("simulations/arrivaldata_etas.rda")

# if (!file.exists("eta_quantiles_A.rda")) {
#     eta_quantiles <-
#         read_csv("simulations/sim002/eta_quantiles.csv",
#             col_names = c("trip_id", "vehicle_id", "stop_sequence", "timestamp", "eta", "quantile"),
#             col_types = "cciiin"
#         )
#     save(eta_quantiles, file = "eta_quantiles_A.rda")
# } else {
#     load("eta_quantiles_A.rda")
# }

## Write CSV to SQLite db
library(RSQLite)
library(dbplyr)
db <- "simulations/sim002/eta_quantiles.sqlite"
con <- dbConnect(SQLite(), db)

if (!dbExistsTable(con, "quantiles")) {
    dbCreateTable(con, "quantiles",
        c(
            trip_id = "text",
            vehicle_id = "text",
            stop_sequence = "int",
            timestamp = "int",
            eta = "int",
            quantile = "real"
        )
    )
    dbExecute(con, "CREATE INDEX tid ON quantiles(trip_id);")
    cat(".mode csv\n.import simulations/sim002/eta_quantiles.csv quantiles\n",
        file = "simulations/sim002/imp.sql")
    system("sqlite3 simulations/sim002/eta_quantiles.sqlite < simulations/sim002/imp.sql")
    unlink("simulations/sim002/imp.sql")
}

# clean up:
# dbRemoveTable(con, "etas")
if (!dbExistsTable(con, "etas")) {
    dbCreateTable(con, "etas",
        c(
            trip_id = "text",
            vehicle_id = "text",
            stop_sequence = "int",
            timestamp = "int",
            eta = "int",
            quantile = "real",
            scheduled_arrival = "int",
            actual_arrival = "int",
            delay = "int",
            time_until_arrival = "int",
            current_delay = "int",
            gtfs_arrival_time = "int",
            gtfs_eta = "int"
        )
    )
    dbExecute(con, "CREATE INDEX tid ON etas(trip_id);")
}

trips <- con %>% tbl("quantiles") %>% select(trip_id) %>% distinct() %>% collect() %>% pull("trip_id")
trips_complete <- con %>% tbl("etas") %>% select(trip_id) %>% distinct() %>% collect() %>% pull("trip_id")
trips <- trips[!trips %in% trips_complete]
pbapply::pblapply(trips, function(trip) {
    eta_quantiles <- con %>%
        tbl("quantiles") %>%
        filter(trip_id == !!trip) %>%
        collect() %>%
        left_join(
            arrivaldata %>%
                filter(trip_id == trip) %>%
                rename(actual_arrival = arrival_time) %>%
                select(trip_id, vehicle_id, stop_sequence, scheduled_arrival, actual_arrival, delay),
            by = c("trip_id", "vehicle_id", "stop_sequence")
        ) %>%
        mutate(
            scheduled_arrival = as.integer(scheduled_arrival),
            actual_arrival = as.integer(actual_arrival),
            time_until_arrival = actual_arrival - timestamp
        ) %>%
        group_by(vehicle_id) %>%
        group_modify(~ {
            # .y <- tibble(trip_id = unique(etas$trip_id)[2])
            # .x <- eta_data %>% filter(trip_id == .y$trip_id)

            d <- rawarrival %>%
                filter(
                    trip_id == trip &
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
        filter(time_until_arrival >= 0)
    dbWriteTable(con, "etas", eta_quantiles, append = TRUE)
})

eta_quantiles <- con %>% tbl("etas") %>%
    select(trip_id, vehicle_id, stop_sequence, timestamp,
        eta, quantile, actual_arrival, time_until_arrival, gtfs_eta)


# eta_quantiles <- eta_quantiles %>%
#     mutate(timestamp = ts2dt(timestamp))
trips <- con %>% tbl("etas") %>% select(trip_id) %>% distinct() %>% collect() %>% pull(trip_id)

smry_file <- "simulations/sim002/eta_smry.rda"
if (file.exists(smry_file)) {
    load(smry_file)
} else {
    eta_smry <- pbapply::pblapply(trips,
        function(trip) {
            eta_quantiles %>% filter(trip_id == !!trip) %>%
                collect() %>%
                mutate(timestamp = ts2dt(timestamp)) %>%
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
        }
    ) %>% bind_rows %>% ungroup()
    save(eta_smry, file = smry_file)
}

ggplot(eta_smry, aes(stop_sequence, err50)) + geom_point()

i <- 5000
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
    geom_vline(xintercept = floor(eta1$gtfs_eta[1] / 60), colour = "magenta")
