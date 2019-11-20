library(tidyverse)

nw_state <- read_csv("simulations/sim000/segment_states.csv",
    col_names = c("segment_id", "timestamp", "avg_speed", "uncertainty"),
    col_types = "iinn"
) %>% mutate(timestamp = as.POSIXct(timestamp, origin = "1970-01-01"))
nw_obs <- read_csv("simulations/sim000/segment_observations.csv",
    col_names = c("segment_id", "timestamp", "speed", "error"),
    col_types = "iinn"
) %>% mutate(timestamp = as.POSIXct(timestamp, origin = "1970-01-01"))

sids <- nw_state %>% pull(segment_id) %>% table %>%
    sort(decreasing = TRUE) %>% names %>% head(20)

sids <- c(38, 46, 97, 161, 187, 155)

ggplot(nw_state %>% filter(segment_id %in% sids), aes(timestamp)) +
    geom_pointrange(
        aes(
            y = speed,
            ymin = truncnorm::qtruncnorm(0.025, 0, 110/3.6, speed, sqrt(error)),
            ymax = truncnorm::qtruncnorm(0.975, 0, 110/3.6, speed, sqrt(error))
        ),
        data = nw_obs %>% filter(segment_id %in% sids)
    ) +
    geom_ribbon(
        aes(
            ymin = truncnorm::qtruncnorm(0.025, 0, 110/3.6, avg_speed, sqrt(uncertainty) + 2),
            ymax = truncnorm::qtruncnorm(0.975, 0, 110/3.6, avg_speed, sqrt(uncertainty) + 2)
        ),
        fill = "gray",
        alpha = 0.5
    ) +
    geom_ribbon(
        aes(
            ymin = truncnorm::qtruncnorm(0.025, 0, 110/3.6, avg_speed, sqrt(uncertainty)),
            ymax = truncnorm::qtruncnorm(0.975, 0, 110/3.6, avg_speed, sqrt(uncertainty))
        ),
        fill = "orangered",
        alpha = 0.5
    ) +
    geom_path(aes(y = avg_speed),
        colour = "orangered"
    ) +
    facet_wrap(~segment_id, ncol = 2) +
    scale_y_continuous(
        name = "speed (m/s)",
        limits = c(0, 110 / 3.6),
        sec.axis = sec_axis(
            ~ . * 3.6,
            name = "speed (km/h)"
        )
    ) -> ppp
dev.hold(); print(ppp); dev.flush(dev.flush())
