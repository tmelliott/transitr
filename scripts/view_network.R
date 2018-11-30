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
    data <- get_segments(f) %>% filter(timestamp > 0)
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


### raw "data"

read_segment_data <- function(sim) {
    file.path("simulations", sim, "history") %>%
        list.files(pattern = "segment_", full.names = TRUE) %>%
        lapply(read_csv, col_types = "ciid", col_names = c("segment_id", "timestamp", "travel_time", "error")) %>%
        bind_rows() %>%
        mutate(timestamp = as.POSIXct(timestamp, origin = "1970-01-01"))
}

## segment lengths?
con <- dbConnect(SQLite(), "fulldata.db")
seglens <- con %>% tbl("road_segments") %>% select(road_segment_id, length) %>% collect %>%
    mutate(road_segment_id = as.character(road_segment_id))
dbDisconnect(con)

segdata <- read_segment_data("sim000")
# segdata <- read_segment_data("sim002")
segids <- table(segdata$segment_id) %>% sort %>% tail(100) %>% names %>% sample(20)

# segids <- table(segdata$segment_id) %>% names %>% sample(20)
segd <- segdata %>% filter(segment_id %in% segids) %>% 
    left_join(seglens, by = c("segment_id" = "road_segment_id")) %>%
    mutate(speed = length / travel_time)

# ggplot(segd) +
#     # geom_smooth(aes(timestamp, travel_time), formula = y ~ 1, method = "lm") +
#     # geom_smooth(aes(timestamp, travel_time)) +
#     geom_hline(aes(yintercept = length / (100 * 1000 / 60 / 60)), color = "orangered") +
#     geom_hline(aes(yintercept = length / (50 * 1000 / 60 / 60)), color = "orangered", lty = 2) +
#     geom_hline(aes(yintercept = length / (30 * 1000 / 60 / 60)), color = "orangered", lty = 3) +
#     geom_hline(aes(yintercept = length / (10 * 1000 / 60 / 60)), color = "orangered", lty = 4) +
#     geom_pointrange(
#         aes(timestamp, travel_time, ymin = pmax(0, travel_time - error), ymax = travel_time + pmin(error, 30)),
#         size = 0.2
#         ) +
#     facet_wrap(~segment_id, scales = "free_y")

ggplot(segd, aes(timestamp, speed / 1000 * 60 * 60)) +
    # geom_smooth(formula = y ~ 1, method = "lm") +
    # geom_smooth(aes(timestamp, travel_time)) +
    geom_hline(aes(yintercept = 10 / 1000 * 60 * 60), color = "orangered") +
    geom_pointrange(
        aes(ymin = pmax(0, length / (travel_time + error) / 1000 * 60 * 60), 
            ymax = pmin(100, length / pmax(1, travel_time - error) / 1000 * 60 * 60)),
        size = 0.2
        ) +
    facet_wrap(~paste0(segment_id, " [", round(length), "m]"))













###############################
library(tidyverse)
library(RSQLite)
library(dbplyr)


shape <- getshape("133")
ggplot(shape[1:100,], aes(shape_pt_lon, shape_pt_lat)) +
    geom_path() +
    geom_point(shape = 4, data = shape %>% filter(shape_pt_sequence == 1)) +
    geom_point(aes(174.635, -36.8782), data = tibble(), color = 'gray') +
    geom_point(aes(174.635, -36.8797), data = tibble()) +
    geom_point(aes(174.636, -36.88), data = NULL, color = "red") +
    geom_point(aes(174.635, -36.8796), data = NULL, color = "orangered")
    
    # geom_point(aes(174.664, -36.8621), data = NULL) +
    # geom_point(aes(174.668, -36.8669), data = NULL, color = "blue") +
    # geom_point(aes(174.646, -36.8706), data = NULL, color = "red")




