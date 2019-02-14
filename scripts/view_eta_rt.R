source("scripts/common.R")

load("simulations/arrivaldata.rda")
arrivaldata <- arrivaldata %>% ungroup()

res <- all_sims("sim000") %>%
    ungroup() %>%
    mutate(dt = as.integer(time - timestamp)) %>%
    filter(dt < 5*60*60) %>%
    mutate(
        lower = as.POSIXct(q0.025, origin = "1970-01-01"),
        upper = as.POSIXct(q0.975, origin = "1970-01-01")
    ) %>%
    group_by(trip_id) %>%
    do((.) %>% mutate(
        tx = as.integer(timestamp),
        tx = (tx - min(tx)) / diff(range(tx))
    )) %>%
    ungroup() %>%
    left_join(
        arrivaldata %>% 
            filter(type == "arrival") %>%
            select(trip_id, stop_sequence, time),
        by = c('trip_id', 'stop_sequence'), 
        suffix = c('', '_actual')
    )


ggplot(res,
    aes(time, stop_sequence + tx*0.8, colour = timestamp)) +
    geom_segment(aes(x = lower, xend = upper, yend = stop_sequence + tx*0.8)) +
    geom_point(
        aes(time_actual, y = stop_sequence + 0.8, colour = NULL, group = trip_id),
        data = res %>% select(time_actual, trip_id, stop_sequence) %>% distinct(),
        colour = "red"
    ) %>%
    facet_wrap(~trip_id)
