library(tidyverse)
library(RSQLite)
library(dbplyr)

get_segments <- function(f = "segment_states.csv") {
    read_csv(f, col_names = c("segment_id", "timestamp", "travel_time", "uncertainty"),
        col_types = "iinn")
}

get_segment_data <- function() {
    con <- dbConnect(SQLite(), "fulldata.db")
    segments <- con %>% tbl("road_segments")
    intersections <- con %>% tbl("intersections")
    segments <- segments %>% 
        inner_join(intersections, by = c("int_from" = "intersection_id"), suffix = c("", "_start")) %>%
        inner_join(intersections, by = c("int_to" = "intersection_id"), suffix = c("", "_end")) %>%
        select(road_segment_id, length, intersection_lat, intersection_lon, intersection_lat_end, intersection_lon_end) %>%
        collect
    dbDisconnect(con)
    segments
}

view_segment_states <- function(f = "segment_states.csv", segment, n = 12) {
    data <- get_segments(f)
    if (missing(segment)) {
        segment <- data %>% distinct %>% pluck("segment_id") %>% 
            table %>% sort %>% tail(n) %>% names
    }
    data <- data %>% filter(segment_id %in% segment) %>% distinct %>%
        mutate(timestamp = as.POSIXct(timestamp, origin = "1970-01-01")) %>%
        arrange(timestamp)

    ggplot(data, aes(timestamp, travel_time)) + 
        geom_linerange(aes(ymin = travel_time - uncertainty, ymax = travel_time + uncertainty)) +
        geom_point() +
        xlab("Time") + ylab("Travel Time (seconds)") +
        facet_wrap(~segment_id)
}

map_segments <- function(f = "segment_states.csv", t = max(data$timestamp)) {
    segdata <- get_segment_data()
    data <- get_segments(f) %>% distinct()
    data <- data %>% filter(timestamp <= t) %>%
        group_by(segment_id) %>%
        do((.) %>% filter(timestamp == max(.$timestamp)))
    data <- segdata %>% 
        inner_join(data, by = c("road_segment_id" = "segment_id")) %>%
        mutate(speed = length / travel_time) %>%
        filter(speed < 35)
    
    ggplot(data, aes(intersection_lon, intersection_lat, colour = speed / 1000 * 60 * 60)) +
        geom_segment(aes(xend = intersection_lon_end, yend = intersection_lat_end)) +
        coord_fixed(1.2) +
        labs(colour = "Speed (km/h)") +
        scale_colour_viridis_c(option = "B", limits = c(0, 100)) +
        xlab("") + ylab("")
}

view_segment_states("simulations/sim000/segment_states.csv")

map_segments()
