source("scripts/common.R")


## load vehicle sim data
get_vehicle_data <- function(vid) {
    read_csv(
        file.path("simulations", "sim000", "history", sprintf("v%s_mutate.csv", vid)),
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
