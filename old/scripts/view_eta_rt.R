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


## can we calculate the current PREDICTION ERROR
# i.e., predicted - actual
# for this, we need to know the SCHEDULED ARRIVAL TIMES
get_stop_times <- function(tid) {
    con <- dbConnect(SQLite(), "fulldata.db")
    on.exit(dbDisconnect(con))
    con %>% tbl("stop_times") %>% filter(trip_id == tid) %>%
        select(stop_sequence, arrival_time) %>% arrange(stop_sequence) %>% collect
}
etatbl <- res %>% 
    group_by(trip_id) %>%
    do({
        d <- .$timestamp[1] %>% format("%Y-%m-%d")
        ta <- get_stop_times((.)$trip_id[1]) %>%
            mutate(arrival_time = as.POSIXct(paste(d, arrival_time)))
        (.) %>% left_join(ta, by = "stop_sequence")
    })

## now need to calculate delay as a function of time
etatbl <- etatbl %>% 
    mutate(sched_delay = as.integer(arrival_time - time_actual)) %>%
    filter(sched_delay > -30*60 & sched_delay < 60*60)

etatbl %>% ggplot(aes(sched_delay/60)) + geom_histogram()

x <- etatbl %>% filter(trip_id == etatbl$trip_id[10]) %>%
    filter(time_actual >= timestamp)

    do({
        x <- (.) %>% filter(time_actual >= timestamp)
        st <- bind_rows(
            tibble(stop_sequence = 0, delay = NA, arr = min(x$timestamp)),
            x %>% group_by(stop_sequence) %>%
                summarize(delay = last(sched_delay), arr = last(time_actual))
            )
        x$delay <- st$delay[sapply(x$timestamp, function(t) which(st$arr <= t) %>% max)]
        ggplot(x, aes(time_actual + delay, timestamp)) + 
            geom_point() + 
            geom_point(aes(x = time_actual + as.integer(time - time_actual)), col = 'red') +
            geom_point(aes(x = time_actual), size = 0.5, colour = "blue") +
            facet_wrap(~stop_sequence, scales = "free_x")
        
    })
