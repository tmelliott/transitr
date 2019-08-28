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
    con <- dbConnect(SQLite(), "at_gtfs.db")
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
    nodes <- con %>% tbl("nodes")
    segments <- segments %>% 
        inner_join(nodes, by = c("node_from" = "node_id"), suffix = c("", "_start")) %>%
        inner_join(nodes, by = c("node_to" = "node_id"), suffix = c("", "_end")) %>%
        select(road_segment_id, length, node_lat, node_lon, node_lat_end, node_lon_end) %>%
        collect
    segments
}

view_segment_states <- function(f = "segment_states.csv", o, segment, n = 12, speed = FALSE, id) {
    data <- get_segments(f) %>% filter(timestamp > 0)
    if (!missing(id)) data <- data %>% filter(segment_id %in% id)
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
    # yr <- range(data$.y) * c(1, 1.5)
    # yr <- extendrange(data$.y, f = 0.5)
    p <- ggplot(data, aes(timestamp, .y)) + 
        geom_linerange(aes(ymin = .y - .e, ymax = .y + .e),  color = 'gray') +
        geom_point() +
        xlab("Time") + 
        ylab(ifelse(speed, "Speed (km/h)", "Travel Time (seconds)")) +
        facet_wrap(~segment_id, scales = ifelse(speed, "fixed", "free_y"))
    if (show.data) {
        # yr <- extendrange(c(data$.y, obs$.y))
        p <- p + 
            # geom_linerange(aes(ymin = .y - .e, ymax = .y + .e),  
            #     data = obs,
            #     color = 'pink', lwd = 0.5) +
            geom_point(data = obs, colour = "red", cex = 0.5)
    }
    if (speed) p <- p + ylim(0, 100) #else p <- p + ylim(yr[1], yr[2])
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

    p <- ggplot(data, aes(node_lon, node_lat, 
                     xend = node_lon_end, yend = node_lat_end)) +
        geom_segment(data = segdata, colour = "black", alpha = 0.05) + 
        geom_segment(aes(color = speed)) +
        coord_fixed(1.2) +
        labs(colour = "Speed (km/h)") +
        theme(legend.position = "bottom") +
        scale_colour_viridis_c(option = "B", limits = c(0, 100)) +
        facet_grid(~speed_fct) +
        xlab("") + ylab("")
    p
}


view_segment_states("simulations/sim000/segment_states.csv", "simulations/sim000/segment_observations.csv", n = 20)
view_segment_states("simulations/sim000/segment_states.csv", "simulations/sim000/segment_observations.csv", speed = TRUE, n = 20)
view_segment_states("simulations/sim000/segment_states.csv", n = 20)
view_segment_states("simulations/sim000/segment_states.csv", speed = TRUE, n = 20)

## what is happening in segment 4691
view_segment_states("simulations/sim000/segment_states.csv", "simulations/sim000/segment_observations.csv", n = 20, id = "4691")
library(leaflet)
sx <- get_segment_data() %>% filter(road_segment_id == "4691")

m <- leaflet() %>%
    addTiles() %>%
    addMarkers(sx$intersection_lon, sx$intersection_lat, 
        popup = "Start"
    ) %>%
    addMarkers(sx$intersection_lon_end, sx$intersection_lat_end,
        popup = "End"
    )

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
con <- dbConnect(SQLite(), "at_gtfs.db")
seglens <- con %>% tbl("road_segments") %>% 
    select(road_segment_id, length, node_from, node_to) %>% collect %>%
    mutate(road_segment_id = as.character(road_segment_id))
dbDisconnect(con)

segdata <- read_segment_data("sim000")
# segdata <- read_segment_data("sim002")
segids <- table(segdata$segment_id) %>% sort %>% tail(50) %>% names

# segids <- table(segdata$segment_id) %>% names %>% sample(20)
segd <- segdata %>% filter(segment_id %in% segids) %>% 
    left_join(seglens, by = c("segment_id" = "road_segment_id")) %>%
    mutate(speed = length / travel_time)


segdf <- segd %>% 
    filter(travel_time >= length / 30 & travel_time <= length / 2)

ggplot(segd, aes(timestamp, speed / 1000 * 60 * 60)) +
    # geom_smooth(formula = y ~ 1, method = "lm") +
    # geom_smooth(aes(timestamp, travel_time)) +
    # geom_hline(aes(yintercept = 10 / 1000 * 60 * 60), color = "orangered") +
    geom_point(
        # aes(ymin = pmax(0, length / (travel_time + error) / 1000 * 60 * 60), 
        #     ymax = pmin(100, length / pmax(1, travel_time - error) / 1000 * 60 * 60)),
        size = 0.5
        ) +
    facet_wrap(~paste0(segment_id, " [", round(length), "m]"))

ggplot(segd, aes(speed / 1000 * 60 * 60)) + geom_histogram() + facet_wrap(~segment_id)

ggplot(segdf) +
    # geom_smooth(aes(timestamp, travel_time), formula = y ~ 1, method = "lm") +
    # geom_smooth(aes(timestamp, travel_time)) +
    # geom_hline(aes(yintercept = length / (100 * 1000 / 60 / 60)), color = "orangered") +
    # geom_hline(aes(yintercept = length / (50 * 1000 / 60 / 60)), color = "orangered", lty = 2) +
    # geom_hline(aes(yintercept = length / (30 * 1000 / 60 / 60)), color = "orangered", lty = 3) +
    # geom_hline(aes(yintercept = length / (10 * 1000 / 60 / 60)), color = "orangered", lty = 4) +
    geom_pointrange(
        aes(timestamp, travel_time, ymin = pmax(0, travel_time - error), ymax = travel_time + pmin(error, 30)),
        size = 0.2
        ) +
    facet_wrap(~paste0(segment_id, " [", round(length), "m]"), scales = "free_y")




segdat <- segdata %>% 
    left_join(seglens, by = c("segment_id" = "road_segment_id")) %>%
    filter(travel_time >= length / 30 & travel_time <= length / 2)

segdat %>%
    group_by(segment_id) %>%
    summarize(tt = mean(travel_time), e = sd(travel_time), n = n()) %>%
    filter(n > 1) %>%
    ggplot(aes(log(tt), log(e))) +
        geom_point()

# format the data 
segdat5 <- segdat %>%
    # filter(segment_id %in% unique(segdat$segment_id)[1:50]) %>%
    arrange(segment_id, timestamp) %>% 
    select(segment_id, timestamp, travel_time, error, length) %>%
    mutate(
        timestamp = as.POSIXct(
            paste0(
                format(timestamp, "%Y-%m-%d %H:%M:"),
                30 * as.integer(format(timestamp, "%S")) %/% 30
            )
        ),
        l = as.integer(as.factor(segment_id)),
        t = as.integer(timestamp)
    ) %>%
    group_by(l) %>%
    do(
        (.) %>% mutate(
            # t0 = as.integer(t == min(t)),
            c = c(1, diff(t) > 0)
        )
    ) %>% ungroup() %>%
    mutate(
        c = cumsum(c)
    )

ggplot(segdat5, aes(timestamp, travel_time)) +
    geom_point() +
    facet_wrap(~segment_id, scales = "free_y")

jdata <-
    list(
        b = segdat5$travel_time,
        e = segdat5$error,
        # identify index of the BETAs
        ell = segdat5$l,
        c = segdat5$c,
        c0 = tapply(segdat5$c, segdat5$l, min) %>% as.integer,
        c1 = tapply(segdat5$c, segdat5$l, min) %>% as.integer + 1,
        cJ = tapply(segdat5$c, segdat5$l, max) %>% as.integer,
        delta = do.call(c, tapply(segdat5$t, segdat5$l, function(tt) {
                c(0, diff(unique(tt)))
            })) %>% as.integer,
        L = length(unique(segdat5$segment_id)),
        N = nrow(segdat5),
        mu = (tapply(segdat5$length, segdat5$l, min) / 30) %>% as.numeric
    )
jdata_t <- do.call(c, tapply(segdat5$timestamp, segdat5$l, unique))
names(jdata_t) <- NULL

model.jags <- "
model{
    for (i in 1:N) {
        b[i] ~ dnorm(B[i], pow(e[i], -2))
        B[i] ~ dnorm(mu[ell[i]] + beta[c[i]], pow(phi[ell[i]], -2))
    }
    
    for (l in 1:L) {
        beta[c0[l]] ~ dunif(0, 500)
        for (j in c1[l]:cJ[l]) {
            beta[j] ~ dnorm(beta[j-1], pow(delta[j] * q[l], -2))
        }
            
        # log(phi) = theta0 + theta1 * log(mu) + err
        phi[l] <- exp(log_phi[l])
        log_phi[l] ~ dnorm(log_mu_phi_ell[l], pow(sig_phi, -2))
        log_mu_phi_ell[l] <- theta[1] + theta[2] * log(mu[l])
    

        q[l] ~ dnorm(mu_q, pow(sig_q, -2))T(0,)
        #log_q[l] ~ dnorm(log_mu_q, pow(sig_q, -2))
        #eq[l] <- exp(log_q[l])
        #q[l] <- pow(eq[l], 2)
    }
    
    theta[1] ~ dnorm(0, 0.01)
    theta[2] ~ dnorm(1, 0.01)
    sig_phi ~ dgamma(0.001, 0.001)
        
    mu_q <- exp(log_mu_q)
    log_mu_q ~ dnorm(0, 0.01)
    sig_q ~ dgamma(0.001, 0.001)
}
"

library(rjags)
library(tidybayes)

jm <- 
    jags.model(textConnection(model.jags),
        # quiet = TRUE
        data = jdata,
        n.chains = 4,
        n.adapt = 10000
    )

n1_samples <-
    coda.samples(jm,
        variable.names = c("beta", "phi", "q", "theta", "sig_phi", "mu_q", "sig_q"),
        n.iter = 5000,
        thin = 5
    )

save(n1_samples, file = "samples_fullmodel.rda")

# n1_samples %>% spread_draws(phi[l]) %>%
#     ggplot() +
#         geom_path(aes(.iteration, phi, colour = as.factor(.chain), group = .chain)) +
#         facet_wrap(~l, scales = "free_y")

n1_samples %>% spread_draws(theta[k]) %>%
    ggplot() +
        geom_path(aes(.iteration, theta, colour = as.factor(.chain), group = .chain)) +
        facet_wrap(~k, scales = "free_y", nrow = 2)


n1_samples %>% spread_draws(q[l]) %>%
    ggplot() +
        geom_path(aes(.iteration, log(q), group = .chain)) +
        facet_wrap(~l, scales = "free_y")

n1_samples %>% spread_draws(mu_q, sig_q) %>%
    ggplot() +
        geom_point(aes(log(mu_q), sig_q, colour = as.factor(.chain))) 


betas <- n1_samples %>% spread_draws(beta[i]) %>%
    median_qi() %>%
    mutate(
        l = tapply(jdata$ell, jdata$c, min),
        t = jdata_t,
        mu = jdata$mu[l]
    )

lii <- table(segdat5$l) %>% sort() %>% tail(20) %>% names()
# lii <- unique(segdat5$l)[1:20 + 80]
ggplot(betas %>% filter(l %in% lii)) +
    geom_ribbon(aes(t, ymin = mu + .lower, ymax = mu + .upper)) +
    geom_path(aes(t, mu + beta)) +
    geom_point(aes(timestamp, travel_time), data = segdat5 %>% filter(l %in% lii), 
        colour = "red", size = 0.5) +
    facet_wrap(~l, scales = "free_y") 

bl <- tapply(jdata$ell, jdata$c, min)
betas_all <- n1_samples %>% spread_draws(beta[i]) %>%
    mutate(l = bl[i], t = jdata_t[i])

ggplot(betas_all) +
    geom_point(aes(t, beta)) +
    geom_point(aes(timestamp, travel_time), data = segdat5, colour = "red") +
    facet_wrap(~l, scales = "free_y")


# segi <- 3316
# ggplot(segdat, aes(intersection_lon, intersection_lat, 
#                  xend = intersection_lon_end, yend = intersection_lat_end)) +
#     geom_segment(colour = "black", alpha = 0.05) + 
#     geom_segment(data = segdat %>% filter(road_segment_id == segi), color = "orangered") +
#     coord_fixed(1.2) +
#     xlab("") + ylab("")


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





### BURP
pf_times <- read_csv("simulations/sim000/particle_travel_times.csv",
    col_types = list(
        col_integer(), col_factor(ordered = TRUE),
        col_double(), col_double()
    ),
    col_names = c("timestamp", "segment_index", "time", "weight")) %>%
    mutate(timestamp = as.POSIXct(timestamp, origin = "1970-01-01"))



### OKAY lets look at coming up with a proper model thing ...
library(rjags)
library(tidybayes)

segdata_all <- get_data("simulations/sim000/segment_observations.csv") %>%
    filter(segment_id %in% 
        (table(segment_id) %>% sort %>% tail(8) %>% names)
    ) %>%
    group_by(segment_id) %>%
    do(
        (.) %>% 
            filter(obs_tt < 5 * median(obs_tt)) %>%
            filter(est_error < 2 * max(obs_tt))
            # mutate(
            #     est_error = pmin(est_error, max(obs_tt))
            # )
    )

ggplot(segdata_all, aes(obs_tt)) + geom_histogram() + facet_wrap(~segment_id)
p0 <- ggplot(segdata_all, aes(timestamp, obs_tt)) +
    geom_pointrange(
        aes(ymin = obs_tt - sqrt(est_error), ymax = obs_tt + sqrt(est_error))
    ) +
    facet_wrap(~segment_id, ncol = 2, scales = "free")
p0

# layout the data for a single segment
segdata <- segdata_all %>% filter(segment_id == 295)
ggplot(segdata, aes(timestamp, obs_tt)) +
    geom_pointrange(
        aes(ymin = obs_tt - sqrt(est_error), ymax = obs_tt + sqrt(est_error))
    )
jags.data <- 
    list(
        # the observations
        b = segdata$obs_tt,
        # the observation errors,
        e = pmax(10, segdata$est_error),
        # the timestamp INDEX
        t = as.integer(as.factor(segdata$timestamp)),
        # the first timestamp (not used yet)
        # T0 = min(segdata$timestamp),
        # time differences
        delta = diff(unique(segdata$timestamp)),
        # number of observations
        N = nrow(segdata),
        # number of timestamps
        M = length(unique(segdata$timestamp))
    )

jags.fit <- 
    jags.model("scripts/nw_model_obs.jags",
        data = jags.data,
        n.chains = 4,
        n.adapt = 5000,
        inits = function() {
            list(
                beta = rnorm(jags.data$M,
                    tapply(segdata$obs_tt, segdata$timestamp, mean),
                    sd(segdata$obs_tt)
                )
                # B = rnorm(jags.data$N, jags.data$b, jags.data$e)
                # log_epsilon = log(runif(1, 1, 5))
            )
        }
    )

samples <- 
    coda.samples(jags.fit,
        # variable.names = c("beta", "kappa", "epsilon", "psi"),
        # variable.names = c("beta", "kappa", "psi", "B"),
        variable.names = c("beta", "B", "kappa", "psi", "epsilon"),
        n.iter = 10000,
        thin = 10
    )



## just the BETA values
beta.samples <- samples %>% spread_draws(beta[t]) %>%
    mutate(timestamp = unique(segdata$timestamp)[t])

ggplot(beta.samples, aes(timestamp, beta)) +
    geom_point(size = 0.2) +
    facet_wrap(~.chain) + 
    geom_pointrange(
        aes(
            y = obs_tt, 
            ymin = obs_tt - sqrt(est_error), 
            ymax = obs_tt + sqrt(est_error)
        ),
        data = segdata,
        colour = "red",
        size = 0.2
    )



samples %>% spread_draws(epsilon) %>% 
    ggplot(aes(.iteration, epsilon, group = .chain, colour = as.factor(.chain))) +
    geom_path() 

samples %>% spread_draws(kappa) %>% 
    ggplot(aes(.iteration, kappa, colour = as.factor(.chain))) + geom_path()

samples %>% spread_draws(kappa[t]) %>% 
    mutate(timestamp = unique(segdata$timestamp)[t]) %>%
    ggplot(aes(timestamp, kappa)) + 
    geom_point() +
    facet_wrap(~.chain)

samples %>% spread_draws(beta[t], kappa[t]) %>% 
    mutate(timestamp = unique(segdata$timestamp)[t]) %>%
    mutate(beta_min = beta - 2*kappa, beta_max = beta + 2*kappa) %>%
    ggplot(aes(timestamp, beta)) +
    geom_linerange(aes(y=NULL, ymin = beta_min, ymax = beta_max)) +
    facet_wrap(~.chain)

     + 
    geom_pointrange(
        aes(
            y = obs_tt, 
            ymin = obs_tt - sqrt(est_error), 
            ymax = obs_tt + sqrt(est_error)
        ),
        data = segdata,
        colour = "red"
    )


samples %>% spread_draws(psi) %>% 
    ggplot(aes(.iteration, psi, group = .chain, colour = as.factor(.chain))) +
    geom_path() 

samples %>% spread_draws(B[i]) %>% 
    mutate(timestamp = segdata$timestamp[i]) %>%
    ggplot(aes(timestamp, B)) +
    geom_point() +
    facet_wrap(~.chain) +
    geom_pointrange(
        aes(
            y = obs_tt, 
            ymin = obs_tt - sqrt(est_error), 
            ymax = obs_tt + sqrt(est_error)
        ),
        data = segdata,
        colour = "red",
        size = 0.5
    )


## Now another model, which uses max speed ...
con <- dbConnect(SQLite(), "at_gtfs.db")
seglens <- con %>% tbl("road_segments") %>% 
    filter(road_segment_id %in% !!unique(segdata$segment_id)) %>%
    collect() %>%
    pluck("length")

jags.data <- 
    list(
        # the observations
        b = segdata$obs_tt,
        # the observation errors,
        e = pmax(10, segdata$est_error),
        # the timestamp INDEX
        t = as.integer(as.factor(segdata$timestamp)),
        # the first timestamp (not used yet)
        # T0 = min(segdata$timestamp),
        # time differences
        delta = diff(unique(segdata$timestamp)),
        # number of observations
        N = nrow(segdata),
        # number of timestamps
        M = length(unique(segdata$timestamp)),
        # segment length
        length = seglens
    )

jags.fit <- 
    jags.model("scripts/nw_model_obs_length.jags",
        data = jags.data,
        inits = function() {
            list(
                #Bmin = runif(1, 40, 60),
                speed_i = 3,
                beta = truncnorm::rtruncnorm(jags.data$M,
                    a = 0, b = Inf,
                    tapply(segdata$obs_tt, segdata$timestamp, mean),
                    sd(segdata$obs_tt)
                )
                # B = rnorm(jags.data$N, jags.data$b, jags.data$e)
                # log_epsilon = log(runif(1, 1, 5))
            )
        },
        n.chains = 4,
        n.adapt = 1000,
    )

samples <- 
    coda.samples(jags.fit,
        # variable.names = c("beta", "kappa", "epsilon", "psi"),
        # variable.names = c("beta", "kappa", "psi", "B"),
        variable.names = c("beta", "B", "kappa", "psi", "epsilon", "Bmin", "max_speed"),
        n.iter = 1000,
        thin = 1
    )

samples %>% spread_draws(Bmin) %>%
    ggplot(aes(.iteration, Bmin, group = .chain, colour = as.factor(.chain))) +
    geom_path() 

samples %>% spread_draws(max_speed) %>%
    ggplot(aes(.iteration, max_speed, group = .chain, colour = as.factor(.chain))) +
    geom_path() 

ggplot(segdata, aes(timestamp, jags.data$length / obs_tt / 1000 * 60 * 60)) +
    geom_point()


## just the BETA values
beta.samples <- samples %>% spread_draws(beta[t], Bmin) %>%
    mutate(timestamp = unique(segdata$timestamp)[t])

ggplot(beta.samples, aes(timestamp, Bmin + beta)) +
    geom_point(size = 0.2) +
    facet_wrap(~.chain) + 
    geom_pointrange(
        aes(
            y = obs_tt, 
            ymin = obs_tt - sqrt(est_error), 
            ymax = obs_tt + sqrt(est_error)
        ),
        data = segdata,
        colour = "red",
        size = 0.2
    )



samples %>% spread_draws(epsilon) %>% 
    ggplot(aes(.iteration, epsilon, group = .chain, colour = as.factor(.chain))) +
    geom_path() 

samples %>% spread_draws(kappa) %>% 
    ggplot(aes(.iteration, kappa, colour = as.factor(.chain))) + geom_path()

samples %>% spread_draws(kappa[t]) %>% 
    mutate(timestamp = unique(segdata$timestamp)[t]) %>%
    ggplot(aes(timestamp, kappa)) + 
    geom_point() +
    facet_wrap(~.chain)

samples %>% spread_draws(beta[t], kappa[t]) %>% 
    mutate(timestamp = unique(segdata$timestamp)[t]) %>%
    mutate(beta_min = beta - 2*kappa, beta_max = beta + 2*kappa) %>%
    ggplot(aes(timestamp, beta)) +
    geom_linerange(aes(y=NULL, ymin = beta_min, ymax = beta_max)) +
    facet_wrap(~.chain)

     + 
    geom_pointrange(
        aes(
            y = obs_tt, 
            ymin = obs_tt - sqrt(est_error), 
            ymax = obs_tt + sqrt(est_error)
        ),
        data = segdata,
        colour = "red"
    )


samples %>% spread_draws(psi) %>% 
    ggplot(aes(.iteration, psi, group = .chain, colour = as.factor(.chain))) +
    geom_path() 

samples %>% spread_draws(B[i]) %>% 
    mutate(timestamp = segdata$timestamp[i]) %>%
    ggplot(aes(timestamp, B)) +
    facet_wrap(~.chain) +
    geom_pointrange(
        aes(
            y = obs_tt, 
            ymin = obs_tt - sqrt(est_error), 
            ymax = obs_tt + sqrt(est_error)
        ),
        data = segdata,
        colour = "red",
        size = 0.5
    ) +
    geom_point(size = 0.2) 

