# estimate delay at start of route
source("scripts/common.R")

load("simulations/arrivaldata.rda")
load("tripupdates.rda")

dep0 <- tripupdates %>% ungroup %>% 
    mutate(date = as.Date(format(timestamp, "%Y-%m-%d"))) %>%
    group_by(vehicle_id, date, trip_id) %>%
    do({
        x <- (.) %>% arrange(time) %>% mutate(dt = c(0, diff(stop_sequence)), n = n())
        # if (all(x$n < 2)) return(NULL)
        if (any(x$dt < 0)) {
            ## need to choose the longest series /only/
            x$ct <- cumsum(x$dt < 0)
            x <- x %>% filter(x$ct == names(sort(table(x$ct), decreasing = TRUE))[1])
            x <- x %>% select(-ct)
        }
        x %>% select(-dt)
    }) %>%
    ungroup() %>% filter(stop_sequence == 1 & dt >= 0) %>%
    select(-stop_sequence, -type)

con <- dbConnect(SQLite(), "fulldata.db")
tids <- unique(dep0$trip_id)
sched0 <- con %>% tbl("stop_times") %>% 
    filter(stop_sequence == 1 & trip_id %in% tids) %>% 
    select(trip_id, departure_time) %>% collect
dbDisconnect(con)

daysahead <- function(time) sapply(strsplit(time, ":"), function(x) as.numeric(x[1]) %/% 24)
convtime <- function(time) {
    sapply(strsplit(time, ":"), function(x) {
        x <- as.numeric(x)
        x[1] <- x[1] %% 24
        sprintf("%02d:%02d:%02d", x[1], x[2], x[3])
    })
}

deps <- dep0 %>%
    filter(n > 2) %>% 
    left_join(sched0, by = "trip_id") %>%
    filter(!is.na(departure_time)) %>%
    mutate(
        time = as.POSIXct(time, origin = "1970-01-01"),
        departure_time = as.POSIXct(paste(date, convtime(departure_time), "NZDT")),
        delay = as.integer(time - departure_time)
    ) %>%
    filter(delay > -60*60 & delay < 60*60*2)

ggplot(deps, aes(delay/60)) + geom_histogram() + xlim(-1, 5) + xlab("Departure delay (min)")

ggplot(deps, aes(time, delay/60/60)) +
    geom_point()


con <- dbConnect(SQLite(), "fulldata.db")
rids <- unique(deproute$route_id)
routeinfo <- con %>% tbl("routes") %>% 
    filter(route_id %in% rids) %>% 
    select(route_id, route_short_name) %>% collect
dbDisconnect(con)

deproute <- deps %>%
    group_by(route_id) %>%
    summarize(
        n = n(),
        delay.med = median(delay, na.rm = TRUE),
        delay.q25 = quantile(delay, 0.25),
        delay.q75 = quantile(delay, 0.75),
        delay.err = sd(delay, na.rm = TRUE) / sqrt(n),
        delay = mean(delay, na.rm = TRUE),
    ) %>%
    mutate(route = as.numeric(as.factor(route_id))) %>%
    left_join(routeinfo, by = "route_id")

depx <- deps %>%
    mutate(route = as.numeric(as.factor(route_id))) %>%
    left_join(routeinfo, by = "route_id")

ggplot(depx, aes(delay / 60, route)) +
    geom_point(size = 0.5) +
    geom_vline(aes(xintercept = 5), col = "red", lty = 2) +
    geom_vline(aes(xintercept = -1), col = "red", lty = 2) +
    xlab("Delay (min)") +
    ggtitle(sprintf(
        "Median delays (with 25%% - 75%% quantiles): %.0f%% \"on time\"",
        mean(deps$delay > -60 & deps$delay < 5*60) * 100
    ))

ggplot(deproute, 
    aes(delay.med/60, route)) + 
    geom_vline(aes(xintercept = 5), col = "red", lty = 2) +
    geom_vline(aes(xintercept = -1), col = "red", lty = 2) +
    geom_segment(aes((delay.q25/60), xend = (delay.q75/60), yend = route)) +
    geom_point() +
    xlab("Delay (min)") +
    geom_label(
        aes(x = delay.q75/60, label = route_short_name),
        data = deproute %>% filter(delay.q75 > 5*60),
        nudge_x = 1, size = 1.8, label.padding = unit(0.1, "lines")
    ) +
    geom_label(
        aes(x = delay.q25/60, label = route_short_name),
        data = deproute %>% filter(delay.q25 < -90),
        nudge_x = -1, size = 1.8, label.padding = unit(0.1, "lines")
    ) +
    ggtitle(sprintf(
        "Median delays (with 25%% - 75%% quantiles): %.0f%% \"on time\"",
        mean(deps$delay > -60 & deps$delay < 5*60) * 100
    ))
