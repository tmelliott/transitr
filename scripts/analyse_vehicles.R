library(tidyverse)
library(RSQLite)
library(dbplyr)
library(patchwork)

vehicles <- list.files("simulations/sim000/history",
    pattern = "vehicle_.{4}.csv",
    full.names = TRUE
)
vi <- vehicles[grepl("5636", vehicles)]
vdata <- read_csv(vi) %>% filter(trip_id == t)

con <- dbConnect(SQLite(), "at_gtfs.db")
route_shapes <- con %>% tbl("trips") %>%
    left_join(
        con %>% tbl("shapes"),
        by = "shape_id"
    ) %>%
    filter(trip_id %in% !!unique(vdata$trip_id)) %>%
    arrange(shape_pt_sequence) %>%
    select(trip_id, shape_id, shape_pt_lat, shape_pt_lon, shape_pt_sequence) %>%
    collect()
route_stops <- con %>% tbl("stop_times") %>%
    left_join(con %>% tbl("stops"), by = "stop_id") %>%
    filter(trip_id %in% !!unique(vdata$trip_id)) %>%
    arrange(arrival_time) %>%
    select(trip_id, stop_sequence, arrival_time, departure_time, stop_lon, stop_lat, shape_dist_traveled) %>%
    collect()
dbDisconnect(con)


pmap <- ggplot(vdata %>% filter(state_type != "mutate")) +
    # facet_wrap(~trip_id) +
    geom_path(aes(shape_pt_lon, shape_pt_lat),
        data = route_shapes,
        colour = "magenta"
    ) +
    geom_point(aes(stop_lon, stop_lat),
        data = route_stops,
        shape = 21, fill = "white",
        colour = "magenta"
    ) +
    geom_point(aes(event_longitude, event_latitude), shape = 4) +
    # geom_path(aes(event_longitude, event_latitude)) +
    geom_point(aes(vehicle_longitude, vehicle_latitude, colour = event_type),
        size = 1, shape = 3
    ) +
    coord_fixed(ratio = cos(mean(vdata$event_latitude * pi / 180., na.rm = TRUE)))

pdist <- ggplot(vdata %>% filter(state_type != "mutate")) +
    geom_point(aes(event_timestamp, vehicle_distance, colour = event_type))

pmap + pdist




## -- old stuff:
## load vehicle sim data
get_vehicle_data <- function(vid) {
    read_csv(
        file.path("simulations", "sim000", "history", sprintf("vehicle_%s.csv", vid)),
        col_names = c(
            "timestamp", "trip_id", "event_type", "event_lat", "event_lon",
            "event_stop_index", "dist_to_route", "delta", "distance", "velocity",
            "lat", "lon", "dist_between", "sumlh"
        ),
        col_types = "iccnnininnnnnn"
    ) %>%
    left_join(
        read_csv(
            file.path("simulations", "sim000", "history", sprintf("v%s_update.csv", vid)),
            col_names = c(
                "timestamp", "trip_id", "event_type", "distance", "velocity",
                "lat", "lon", "dist_between", "Neff", "resample"
            ),
            col_types = "iccnnnnnni"
        ),
        by = c("timestamp", "trip_id", "event_type"),
        suffix = c("_prior", "_posterior")
    ) %>%
    mutate(
        timestamp = as.POSIXct(timestamp, origin = "1970-01-01"),
        event_type = factor(event_type, levels = c("initialize", "gps", "arrival", "departure"))
    )
}

vehicles <- list.files(file.path("simulations", "sim000", "history"), pattern = "v*_mutate.csv")
vehicles <- gsub("^v|\\_mutate.csv", "", vehicles)
# vehicles <- list.files("simulations/sim000/history", pattern = "vehicle_.{4}.csv")
# vids <- gsub("vehicle_|.csv", "", vehicles)

draw <- function() {
    vdata <- get_vehicle_data(vid <- sample(vehicles, 1)) %>%
        filter(!is.na(event_lat)) %>% arrange(timestamp)
    ggplot(vdata, aes(event_lon, event_lat, group = trip_id)) +
        geom_path() +
        coord_fixed(ratio = 1.2) +
        ggtitle(sprintf("Vehicle %s from %s - %s (%s obs)",
            vid, format(min(vdata$timestamp), "%T"), format(max(vdata$timestamp), "%T"), nrow(vdata)))
}
draw()

## compute dwell times


ggplot(vdata, aes(timestamp, event_type)) +
    geom_point()

ggplot(vdata, aes(dist_to_route, dist_between_posterior)) +
    geom_point() +
    geom_abline()

ggplot(vdata, aes(distance_posterior/1000, velocity_posterior/1000*60*60)) +
    geom_point() +
    xlab("Distance (km)") + ylab("Speed (km/h)")

ggplot(vdata %>% filter(event_type == "gps"), aes(timestamp, distance_prior)) +
    geom_path() +
    geom_path(aes(y = distance_posterior), color = "orangered")

ggplot(vdata %>% filter(event_type == "gps"), aes(timestamp, dist_to_route)) +
    geom_path() +
    geom_point(aes(y = dist_between_posterior), color = "orangered")

shape <- getshape(trip = vdata$trip_id[1])

vdatagps <- vdata %>% filter(event_type == "gps")

for (t in vdatagps$timestamp %>% unique) {
    p <- ggplot(shape, aes(shape_pt_lon, shape_pt_lat)) +
        geom_path() +
        geom_point(aes(event_lon, event_lat), data = vdatagps %>% filter(timestamp == t)) +
        # geom_point(aes(lon_prior, lat_prior), data = vdatagps, color = "orangered") +
        geom_point(aes(lon_posterior, lat_posterior), data = vdatagps %>% filter(timestamp == t),
            color = "orangered")
    print(p)
    grid::grid.locator()
}
