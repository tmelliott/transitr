source("scripts/common.R")

sim <- "sim000"
vids <- unique(gsub(
    "^v|_mutate\\.csv", 
    "", 
    list.files(
        file.path("simulations", sim, "history"),
        pattern = "^v[A-Z0-9]*_mutate\\.csv"
    )
))

v_hist <- function(id, type = c("mutate", "update"), sim = "sim000") {
    type <- match.arg(type)
    f <- sprintf(
        "simulations/%s/history/v%s_%s.csv",
        sim, id, type
    )
    if (type == "mutate") {
        cn <- c("timestamp", "trip_id", "type", "obs_latitude", "obs_longitude",
            "stop_index", "dist_to_route", "delta", "distance", "speed",
            "latitude", "longitude", "dist_to_obs", "sum_llh")
        ct <- "iccnnininnnnnn"
    } else {
        cn <- c("timestamp", "trip_id", "type", "distance", "speed", "latitude", 
            "longitude", "dist_to_obs", "Neff", "keep")
        ct <- "iccnnnnnni"
    }
    read_csv(f, col_names = cn, col_types = ct)
}

vi <- 2
vv <- v_hist(vids[vi]) %>% 
    left_join(
        v_hist(vids[vi], "update") %>% 
            select(timestamp, trip_id, type, distance, speed, latitude, 
                longitude, dist_to_obs, Neff, keep),
        by = c("timestamp", "trip_id", "type"),
        suffix = c("", "_update")
    ) %>%
    mutate(time = as.POSIXct(timestamp, origin = "1970-01-01"))


ggplot(vv, aes(time)) +
    geom_point(aes(y = distance_update, colour = type, size = sum_llh)) +
    facet_wrap(~trip_id, scales = "free")

ggplot(vv, aes(time)) +
    geom_point(aes(y = dist_to_obs_update, colour = type, size = sum_llh)) +
    facet_wrap(~trip_id, scales = "free")

ggplot(vv %>% filter(type == "gps"), aes(time)) +
    geom_point(aes(y = dist_to_route)) +
    geom_point(aes(y = dist_to_obs_update, size = sum_llh), 
        colour = "orangered", alpha = 0.5) +
    facet_wrap(~trip_id, scales = "free")

sh <- do.call(
    bind_rows,
    lapply(unique(vv$trip_id), function(t) 
        getshape(trip = t) %>% mutate(trip_id = t))
)

ggplot(vv, aes(obs_longitude, obs_latitude)) +
    geom_path(aes(shape_pt_lon, shape_pt_lat, group = shape_id), data = sh) +
    geom_point() +
    geom_segment(aes(xend = longitude, yend = latitude), colour = "orangered") +
    geom_point(aes(longitude, latitude, colour = "orangered"), size = 1) +
    geom_segment(aes(xend = longitude_update, yend = latitude_update), colour = "blue") +
    geom_point(aes(longitude_update, latitude_update), colour = "blue", size = 1) +
    facet_wrap(~trip_id)




vi <- 2
vv2 <- v_hist(vids[vi], sim = "sim001") %>% 
    full_join(
        v_hist(vids[vi], "update", sim = "sim001") %>% 
            select(timestamp, trip_id, type, distance, speed, latitude, 
                longitude, dist_to_obs, Neff, keep),
        by = c("timestamp", "trip_id", "type"),
        suffix = c("", "_update")
    ) %>%
    mutate(time = as.POSIXct(timestamp, origin = "1970-01-01"))

vv2 %>% filter(type == "gps") %>% pluck("keep")

ggplot(vv2, aes(time)) +
    geom_point(aes(y = distance_update, colour = type, size = sum_llh)) +
    facet_wrap(~trip_id, scales = "free")

ggplot(vv2, aes(time)) +
    geom_point(aes(y = dist_to_obs_update, colour = type, size = sum_llh)) +
    facet_wrap(~trip_id, scales = "free")

ggplot(vv2 %>% filter(type == "gps"), aes(time)) +
    geom_point(aes(y = dist_to_route)) +
    geom_point(aes(y = dist_to_obs_update, size = sum_llh), 
        colour = "orangered", alpha = 0.5) +
    facet_wrap(~trip_id, scales = "free")

sh <- do.call(
    bind_rows,
    lapply(unique(vv2$trip_id), function(t) 
        getshape(trip = t) %>% mutate(trip_id = t))
)

ggplot(vv2, aes(obs_longitude, obs_latitude)) +
    geom_path(aes(shape_pt_lon, shape_pt_lat, group = shape_id), data = sh) +
    geom_point() +
    geom_segment(aes(xend = longitude, yend = latitude, alpha = delta), colour = "orangered") +
    geom_point(aes(longitude, latitude, colour = "orangered"), size = 1) +
    geom_segment(aes(xend = longitude_update, yend = latitude_update), colour = "blue") +
    geom_point(aes(longitude_update, latitude_update), colour = "blue", size = 1) +
    facet_wrap(~trip_id)

