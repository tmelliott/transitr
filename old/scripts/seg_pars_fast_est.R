# 'quickly' estimate segment parameters for all segments (mean + se + variance)
library(tidyverse)

load("scripts/data_week0.rda")
load("scripts/data_week1.rda")

tt_data <- bind_rows(data_week0, data_week1) %>%
    mutate(
        time_sec = as.integer(time) - as.integer(as.POSIXct("2019-08-12 00:00:00"))
    )

sid <- unique(tt_data$segment_id)[1]

tt_data %>% filter(segment_id == sid) %>%
    ggplot(aes(time, travel_time)) +
    geom_point(aes(colour = date)) +
    geom_smooth(aes(colour = date)) +
    facet_wrap(~dow)

# compare individual smooths to single linear regression line
library(mgcv)

fit0 <- lm(travel_time ~ dow, data = tt_data)
