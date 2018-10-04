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

view_segment_states <- function(f = "segment_states.csv", segment, n = 12, speed = FALSE) {
    data <- get_segments(f)
    if (missing(segment)) {
        segment <- data %>% distinct %>% pluck("segment_id") %>% 
            table %>% sort %>% tail(n) %>% names
    }
    data <- data %>% filter(segment_id %in% segment) %>% distinct %>%
        mutate(timestamp = as.POSIXct(timestamp, origin = "1970-01-01")) %>%
        arrange(timestamp)

    if (speed) {
        segdata <- get_segment_data() %>% select(road_segment_id, length)
        data <- data %>% 
            inner_join(segdata, by = c("segment_id" = "road_segment_id")) %>%
            mutate(speed = length / travel_time) %>%
            filter(speed < 35) %>%
            mutate(speed = speed * 3.6, .y = speed, .e = sqrt(uncertainty * 3.6^2))
    } else {
        data <- data %>% mutate(.y = travel_time, .e = sqrt(uncertainty))
    }
    p <- ggplot(data, aes(timestamp, .y)) + 
        geom_linerange(aes(ymin = .y - .e, ymax = .y + .e),  color = 'gray') +
        geom_point() +
        xlab("Time") + 
        ylab(ifelse(speed, "Speed (km/h)", "Travel Time (seconds)")) +
        facet_wrap(~segment_id, scales = ifelse(speed, "fixed", "free_y"))
    if (speed) p <- p + ylim(0, 100)
    p
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
        filter(speed < 35) %>%
        mutate(speed = speed / 1000 * 60 * 60,
               speed_fct = case_when(speed < 30 ~ "< 30 kmh",
                                     speed < 55 ~ "30-55 kmh",
                                     speed < 70 ~ "55-70 kmh",
                                     TRUE ~ "70+ kmh"))
    
    ggplot(data, aes(intersection_lon, intersection_lat, 
                     xend = intersection_lon_end, yend = intersection_lat_end)) +
        geom_segment(data = segdata, colour = "black", alpha = 0.05) + 
        geom_segment(aes(color = speed)) +
        coord_fixed(1.2) +
        labs(colour = "Speed (km/h)") +
        theme(legend.position = "bottom") +
        scale_colour_viridis_c(option = "B", limits = c(0, 100)) +
        facet_grid(~speed_fct) +
        xlab("") + ylab("")
}


view_segment_states("simulations/sim000/segment_states.csv", n = 20)
view_segment_states("simulations/sim000/segment_states.csv", speed = TRUE, n = 20)

map_segments("simulations/sim000/segment_states.csv")

view_segment_states("simulations/sim100/segment_states.csv", n = 20)
view_segment_states("simulations/sim100/segment_states.csv", speed = TRUE, n = 20)

map_segments("simulations/sim100/segment_states.csv")
