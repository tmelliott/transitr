library(tidyverse)
library(ggraph)
library(tidygraph)
library(Matrix)
library(rjags)
library(RSQLite)
library(dbplyr)

source("scripts/common2.R")

sim <- "sim000"

## fetch the data
segdata <- get_segment_data()
nwobs <- 
    read_csv(
        file.path("simulations", sim, "segment_observations.csv"),
        col_names = c("segment_id", "timestamp", "travel_time", "measurement_error"),
        col_types = "iiin"
    ) %>% 
    filter(timestamp > 0) %>%
    mutate(
        timestamp = as.POSIXct(timestamp, origin = "1970-01-01")
    ) %>%
    left_join(
        segdata %>% select(road_segment_id, length), 
        by = c("segment_id" = "road_segment_id")
    ) %>%
    mutate(
        speed = length / travel_time,
        speedkm = speed / 1000 * 60 * 60
    )

## look at the data
seg <- nwobs$segment_id %>% table %>% sort %>% names %>% tail(20)

ps <- ggplot(nwobs %>% filter(segment_id %in% seg)) +
    facet_wrap(~segment_id)

ps + geom_histogram(aes(travel_time))
ps + geom_histogram(aes(speedkm))

ps + geom_point(aes(timestamp, travel_time))
ps + geom_point(aes(timestamp, speedkm))

