library(tidyverse)
library(RSQLite)
library(dbplyr)

get_segments <- function(f = "segment_states.csv") {
    read_csv(f, col_names = c("segment_id", "timestamp", "travel_time", "uncertainty"),
        col_types = "iinn")
}

get_data <- function(f = "segment_observations.csv") {
    read_csv(f, col_names = c("segment_id", "timestamp", "obs_tt", "est_error"),
        col_types = "iinn")
}

get_segment_data <- function(routes) {
    con <- dbConnect(SQLite(), "fulldata.db")
    on.exit(dbDisconnect(con))
    if (missing(routes)) {
        segments <- con %>% tbl("road_segments")
    } else {
        rids <- con %>% tbl("routes") %>%
            filter(route_short_name %in% routes &
                (route_long_name %like% "%To City%" |
                 route_long_name %like% "%To Britomart%")) %>%
            collect %>% pluck("route_id") %>% unique
        sids <- con %>% tbl("trips") %>%
            filter(route_id %in% rids) %>%
            collect %>% pluck("shape_id") %>% unique
        segids <- con %>% tbl("shape_segments") %>% 
            filter(shape_id %in% sids) %>%
            collect %>% pluck("road_segment_id") %>% unique
        segments <- con %>% tbl("road_segments") %>%
            filter(road_segment_id %in% segids)
    }
    intersections <- con %>% tbl("intersections")
    segments <- segments %>% 
        inner_join(intersections, by = c("int_from" = "intersection_id"), suffix = c("", "_start")) %>%
        inner_join(intersections, by = c("int_to" = "intersection_id"), suffix = c("", "_end")) %>%
        select(road_segment_id, length, intersection_lat, intersection_lon, intersection_lat_end, intersection_lon_end) %>%
        collect
    segments
}

view_segment_states <- function(f = "segment_states.csv", o, segment, n = 12, speed = FALSE) {
    data <- get_segments(f) %>% filter(timestamp > 0)
    show.data <- !missing(o)
    if (show.data) obs <- get_data(o) %>% filter(timestamp > 0)
    if (missing(segment)) {
        segment <- data %>% distinct %>% pluck("segment_id") %>% 
            table %>% sort %>% tail(n) %>% names
    }
    data <- data %>% filter(segment_id %in% segment) %>% distinct %>%
        mutate(timestamp = as.POSIXct(timestamp, origin = "1970-01-01")) %>%
        arrange(timestamp)
    if (show.data) {
        obs <- obs %>% filter(segment_id %in% segment) %>% distinct %>%
        mutate(timestamp = as.POSIXct(timestamp, origin = "1970-01-01")) %>%
        arrange(timestamp)
    }

    if (speed) {
        segdata <- get_segment_data() %>% select(road_segment_id, length)
        data <- data %>% 
            inner_join(segdata, by = c("segment_id" = "road_segment_id")) %>%
            mutate(speed = length / travel_time) %>%
            filter(speed < 35) %>%
            mutate(speed = speed * 3.6, .y = speed, .e = uncertainty * 3.6^2)
        if (show.data)
            obs <- obs %>%
                inner_join(segdata, by = c("segment_id" = "road_segment_id")) %>%
                mutate(speed = length / obs_tt) %>%
                filter(speed < 35) %>%
                mutate(speed = speed * 3.6, .y = speed, .e = sqrt(est_error * 3.6^2))
    } else {
        data <- data %>% mutate(.y = travel_time, .e = uncertainty)
        if (show.data)
            obs <- obs %>% mutate(.y = obs_tt, .e = sqrt(est_error))
    }
    yr <- extendrange(data$.y, f = 0.5)
    p <- ggplot(data, aes(timestamp, .y)) + 
        geom_linerange(aes(ymin = .y - .e, ymax = .y + .e),  color = 'gray') +
        geom_point() +
        xlab("Time") + 
        ylab(ifelse(speed, "Speed (km/h)", "Travel Time (seconds)")) +
        facet_wrap(~segment_id, scales = ifelse(speed, "fixed", "free_y"))
    if (show.data)
        p <- p + 
            # geom_linerange(aes(ymin = .y - .e, ymax = .y + .e),  
            #     data = obs,
            #     color = 'pink', lwd = 0.5) +
            geom_point(data = obs, colour = "red", cex = 0.5)
    if (speed) p <- p + ylim(0, 100) else p <- p + ylim(yr[1], yr[2])
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


view_segment_states("simulations/sim000/segment_states.csv", "simulations/sim000/segment_observations.csv", n = 20)
view_segment_states("simulations/sim000/segment_states.csv", "simulations/sim000/segment_observations.csv", speed = TRUE, n = 20)
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
seglens <- con %>% tbl("road_segments") %>% 
    select(road_segment_id, length, int_from, int_to) %>% collect %>%
    mutate(road_segment_id = as.character(road_segment_id))
dbDisconnect(con)

segdata <- read_segment_data("sim000")
# segdata <- read_segment_data("sim002")
segids <- table(segdata$segment_id) %>% sort %>% tail(50) %>% names

# segids <- table(segdata$segment_id) %>% names %>% sample(20)
segd <- segdata %>% filter(segment_id %in% segids) %>% 
    left_join(seglens, by = c("segment_id" = "road_segment_id")) %>%
    mutate(speed = length / travel_time)


ggplot(segd, aes(timestamp, speed / 1000 * 60 * 60)) +
    # geom_smooth(formula = y ~ 1, method = "lm") +
    # geom_smooth(aes(timestamp, travel_time)) +
    geom_hline(aes(yintercept = 10 / 1000 * 60 * 60), color = "orangered") +
    geom_point(
        # aes(ymin = pmax(0, length / (travel_time + error) / 1000 * 60 * 60), 
        #     ymax = pmin(100, length / pmax(1, travel_time - error) / 1000 * 60 * 60)),
        size = 0.5
        ) +
    facet_wrap(~paste0(segment_id, " [", round(length), "m]"))

ggplot(segd, aes(speed / 1000 * 60 * 60)) + geom_histogram() + facet_wrap(~segment_id)

ggplot(segd) +
    # geom_smooth(aes(timestamp, travel_time), formula = y ~ 1, method = "lm") +
    # geom_smooth(aes(timestamp, travel_time)) +
    geom_hline(aes(yintercept = length / (100 * 1000 / 60 / 60)), color = "orangered") +
    geom_hline(aes(yintercept = length / (50 * 1000 / 60 / 60)), color = "orangered", lty = 2) +
    geom_hline(aes(yintercept = length / (30 * 1000 / 60 / 60)), color = "orangered", lty = 3) +
    geom_hline(aes(yintercept = length / (10 * 1000 / 60 / 60)), color = "orangered", lty = 4) +
    geom_pointrange(
        aes(timestamp, travel_time, ymin = pmax(0, travel_time - error), ymax = travel_time + pmin(error, 30)),
        size = 0.2
        ) +
    facet_wrap(~paste0(segment_id, " [", round(length), "m]"), scales = "free_y")


segdat <- get_segment_data()
segi <- 3316
ggplot(segdat, aes(intersection_lon, intersection_lat, 
                 xend = intersection_lon_end, yend = intersection_lat_end)) +
    geom_segment(colour = "black", alpha = 0.05) + 
    geom_segment(data = segdat %>% filter(road_segment_id == segi), color = "orangered") +
    coord_fixed(1.2) +
    xlab("") + ylab("")


################ networkisation
library(ggraph)
library(tidygraph)
library(Matrix)

## just segments of these routes
segdata <- get_segment_data(c("27W", "27H", "22N", "22R", "22A"))
sids <- unique(segdata$road_segment_id)

con <- dbConnect(SQLite(), "fulldata.db")
seglens <- con %>% tbl("road_segments") %>% 
    filter(road_segment_id %in% sids) %>%
    select(road_segment_id, length, int_from, int_to) %>% collect %>%
    mutate(road_segment_id = as.character(road_segment_id),
           from = int_from, to = int_to)
dbDisconnect(con)

rawdata <- read_segment_data("sim000") %>% 
    filter(segment_id %in% sids) %>%
    left_join(seglens, by = c("segment_id" = "road_segment_id"))
segg <- rawdata %>%
    mutate(speed = length / travel_time) %>%
    group_by(segment_id) %>% 
        summarize(
            speed = mean(speed) * 60 * 60 / 1000,
            from = first(from), to = first(to)
        ) %>%
    as_tbl_graph()

ggraph(segg, layout = 'kk') + 
    geom_edge_fan(
        aes(color = speed),
        arrow = arrow(length = unit(1, 'mm')), 
        start_cap = circle(1, 'mm'),
        end_cap = circle(1, 'mm')) +
    scale_edge_colour_gradientn(colours = viridis::viridis(256, option = "A"), limits = c(0, 100)) +
    geom_node_point(shape = 19, size = 0.5)

### by the hour (for now)
discretise <- function(x, seconds = 3600) {
    date <- as.integer(as.POSIXct(format(x$timestamp[1], "%Y-%m-%d 00:00:00")))
    bins <- seq(5*60*60, 24*60*60, by = seconds)
    ts <- as.integer(x$timestamp) - date
    x$time <- cut(ts, breaks = bins, 
        labels = as.POSIXct(date + bins[-length(bins)], origin = "1970-01-01") %>%
            format("%T"))
    x
}

seg.hour <- rawdata %>% discretise(60*15) %>%
    mutate(speed = length / travel_time) %>%
    group_by(segment_id, time) %>% 
        summarize(
            speed = mean(speed) * 60 * 60 / 1000,
            from = first(from), to = first(to)
        )
segg.hour <- seg.hour %>% as_tbl_graph() %>% activate(edges)

ggraph(segg.hour, layout = 'kk') + 
    geom_edge_fan(aes(color = speed)) +
    scale_edge_colour_gradientn(colours = viridis::viridis(256, option = "A"), limits = c(0, 100)) +
    facet_wrap(~time)

## smooth them out
segY <- seg.hour %>% spread(key = time, value = speed) %>% ungroup() %>%
    mutate(segment_id = as.character(segment_id)) %>% arrange(segment_id)
Y <- segY %>% select(-from, -to, -segment_id) %>% as.matrix
X <- Y * NA

## should flip Y -> [time x segment]

generate_matrix <- function(x) {
    m <- Matrix(diag(nrow(x)), doDiag = FALSE)
    ## convert id to factor
    x$id <- as.integer(as.factor(x$id))
    ids <- unique(c(x$from, x$to))
    x$from <- as.integer(factor(x$from, levels = ids))
    x$to <- as.integer(factor(x$to, levels = ids))
    for (i in seq_along(1:nrow(m))) {
        toi <- x$id[x$to == x$from[i]]
        fromi <- x$id[x$from == x$to[i]]
        if (length(toi) > 0) m[i, toi] <- 1
        if (length(fromi) > 0) m[i, fromi] <- 1
    }
    sweep(m, 1, rowSums(m), "/")
}

segm <- segY %>% mutate(id = segment_id) %>% select(id, from, to)
F <- Matrix(generate_matrix(segm)) # transition matrix
Q <- Matrix(diag(nrow(X))) * 10
R <- Matrix(diag(nrow(X))) * 5
I <- Matrix(diag(nrow(X)))

xhat <- cbind(rep(30, nrow(X)))
P <- Matrix(diag(nrow(X)) * 100, doDiag = FALSE)

pb <- txtProgressBar(0, ncol(Y), style = 3)
for (i in 1:ncol(Y)) {
    ## predict
    xhat <- as.matrix(F %*% xhat)
    P <- F %*% P %*% t(F) + Q

    ## update
    y <- Y[, i] - xhat
    y[is.na(y)] <- 0
    S <- R + P
    K <- P %*% solve(S)
    X[, i] <- xhat <- xhat + as.matrix(K %*% y)
    P <- (I - K) %*% P %*% t(I - K) + K %*% R %*% t(K)
    setTxtProgressBar(pb, i)
}
close(pb)

segd.orig <- segY %>% gather(key = time, value = speed, -(1:3))
colnames(X) <- names(segY)[-(1:3)]
segd.kf <- segY %>% select(segment_id) %>% bind_cols(as.tibble(X)) %>% 
    gather(key = time, value = speed_kf, -1)

seg.hour <- segd.orig %>% inner_join(segd.kf, by = c('segment_id', 'time'))
segg.hour <- seg.hour %>% as_tbl_graph() %>% activate(edges)

set.seed(1234)
for (t in unique(seg.hour$time)) {
    gt <- ggraph(segg.hour %>% filter(time == t), layout = "kk") + 
        geom_edge_fan(
            aes(color = speed_kf),
            # arrow = arrow(length = unit(1, 'mm')),
            width = 2,
            start_cap = circle(1, 'mm'),
            end_cap = circle(1, 'mm')) +
        scale_edge_colour_gradientn(colours = viridis::viridis(256, option = "A"), limit = c(0, 100)) +
        geom_node_point(shape = 19, size = 1) +
        ggtitle(t)
    dev.hold()
    print(gt)
    dev.flush()
}




animation::saveGIF(
    for(g in gs) {
        dev.hold()
        print(g)
        dev.flush()
    },
    "network.gif",
    ani.width = 800, ani.height = 500
    )


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




