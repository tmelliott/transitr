library(tidyverse)
library(RProtoBuf)

# curd <- setwd("src/vendor/protobuf")
# readProtoFiles("gtfs-realtime.proto")
# setwd(curd)

## load trip updates
if (!file.exists("tripupdates.rda")) {
    tufiles <- list.files(pattern = "trip_updates_.*\\.rda")
    tripupdates <- NULL
    for (f in tufiles) {
        load(f)
        tripupdates <- tripupdates %>% bind_rows(tu %>% distinct())
        rm(tu)
    }
    save(tripupdates, file = "tripupdates.rda")
} else {
    load("tripupdates.rda")
}

## load vehilce positions
# if (!file.exists("vehiclepositions.rda")) {
#     vpfiles <- list.files(file.path("simulations", "archive"), 
#         pattern = "vehicle_locations*", full.name = TRUE)
#     vp <- pbapply::pblapply(vpfiles, function(f) {
#         feed <- read(transit_realtime.FeedMessage, f)
#         if (length(feed$entity) == 0) return(NULL)
#         lapply(feed$entity, function(e) {
#             if (!e$has('vehicle')) return(NULL)
#             if (!e$vehicle$has('position')) return(NULL)
#             xdf <- tibble(
#                 vehicle_id = e$vehicle$vehicle$id,
#                 trip_id = e$vehicle$trip$trip_id,
#                 route_id = e$vehicle$trip$route_id,
#                 timestamp = as.POSIXct(e$vehicle$timestamp, origin = "1970-01-01"),
#                 position_latitude = e$vehicle$position$latitude,
#                 position_longitude = e$vehicle$position$longitude
#             )
#         }) %>% bind_rows
#     }) %>% bind_rows
#     save(vp, file = "vehiclepositions.rda")
# } else {
#     load("vehiclepositions.rda")
# }

dt <- tripupdates %>%
    filter(type == "arrival") %>%
    mutate(arrival = time) %>%
    select(vehicle_id, trip_id, timestamp, stop_sequence, arrival) %>%
    distinct() %>% 
    full_join(
        tripupdates %>% filter(type == "departure") %>% mutate(departure = time) %>%
            select(vehicle_id, trip_id, stop_sequence, departure),
        # by = c("vehicle_id", "trip_id", "stop_sequence"),
        suffix = c(".arrival", ".departure")
    ) %>%
    arrange(vehicle_id, trip_id, stop_sequence) %>%
    mutate(dwell = departure - arrival) %>% filter(!is.na(dwell)) %>%
    filter(dwell >= 0 & dwell <= 1*60)

hp <- ggplot(dt, aes(dwell)) + 
    geom_histogram(binwidth = 1) + 
    xlab("Dwell time (sec)") + ylab("Number of buses") 

hp + geom_vline(aes(xintercept = x), 
    data = tibble(x = c(9, 19, 27, 36, 45)), color = "orangered")

hp + geom_vline(aes(xintercept = x), 
    data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) 


## how about the rate of vehicle updates?
# dv <- vp %>% group_by(vehicle_id) %>% distinct() %>%
#     do((.) %>% arrange(timestamp) %>% mutate(dt = c(-1, diff(timestamp)))) %>%
#     filter(dt >= 0 & dt <= 1*60)

# dp <- ggplot(dv, aes(dt)) +
#     geom_histogram(binwidth = 1) +
    # xlab("Time between updates (sec)") + ylab("Number of updates")

# dp + geom_vline(aes(xintercept = x), 
#     data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) 

# dt2 <- tu %>% group_by(vehicle_id) %>% distinct() %>%
#     do((.) %>% arrange(timestamp) %>% mutate(dt = c(-1, diff(timestamp)))) %>%
#     filter(dt >= 0 & dt <= 1*60)

# hp2 <- ggplot(dt2, aes(dt)) +
#     geom_histogram(binwidth = 1) +
#     xlab("Time between updates (sec)") + ylab("Number of updates")

# hp2 + geom_vline(aes(xintercept = x), 
#     data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) 


## join
# dd <- dv %>% full_join(dt, by = c("vehicle_id", "timestamp")) %>%
#     distinct() %>%
#     do((.) %>% arrange(timestamp) %>% mutate(dt = c(-1, diff(timestamp)))) %>%
#     filter(dt >= 0 & dt <= 1*60)

# dp2 <- ggplot(dd %>% filter(dt < 1*60), aes(dt)) +
#     geom_histogram(binwidth = 1) +
#     xlab("Time between updates (sec)") + ylab("Number of updates")

# dp2 + geom_vline(aes(xintercept = x), 
#     data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) 

# gridExtra::grid.arrange(
#     hp + geom_vline(aes(xintercept = x), 
#         data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) +
#         ggtitle("Dwell times"),
#     hp2 + geom_vline(aes(xintercept = x), 
#         data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) +
#         ggtitle("Trip updates only"),
#     dp + geom_vline(aes(xintercept = x), 
#         data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) +
#         ggtitle("Vehicle positions only"),
#     dp2 + geom_vline(aes(xintercept = x), 
#         data = tibble(x = 9*(1:6)), color = "steelblue", lty = 2) +
#         ggtitle("Combined"),
#     nrow = 4)



# ## find the stop
# library(RSQLite)
# library(dbplyr)
# tids <- unique(dt$trip_id)
# sts <- dbConnect(SQLite(), "fulldata.db") %>% tbl("stop_times") %>% 
#     select(trip_id, stop_sequence, stop_id) %>%
#     filter(trip_id %in% tids) %>%
#     collect

# dt <- dt %>% left_join(sts) %>%
#     mutate(
#         stop_id = as.factor(stop_id),
#         stopid = as.numeric(stop_id)
#     )

# # library(biglm)

# # i <- 1
# # size <- 1000
# # fit <- biglm(dwell ~ stop_id, data = dt[i:(i+size-1),])
# # i <- i+size
# # cat("\n")
# # while (i < nrow(dt)) {
# #     cat(sprintf("\r%i / %i ", i, nrow(dt)))
# #     fit <- update(fit, dt[i:(i+size-1),])
# #     i <- i+size
# # }
# # cat("\r done                  \n")



### Look at travel times

trip_times_to_travel_times <- function(x) {
    ## sort & remove duplicates
    x <- x %>% distinct() %>% arrange(stop_sequence, type) %>%
        mutate(time = as.POSIXct(time, origin = "1970-01-01"))
    ## all should be increasing
    if (nrow(x) < 3 || any(diff(x$time) < 0)) 
        return(x %>% filter(FALSE) %>% 
            mutate(segment_index = stop_sequence, 
                time_from = time, time_to = time, travel_time = numeric()) %>%
            select(route_id, trip_id, segment_index, time_from, time_to, travel_time))
    ## rearrange into from-to
    xFrom <- x %>% 
        filter(stop_sequence < max(stop_sequence) & type == "departure") %>%
        mutate(segment_index = stop_sequence)
    xTo <- x %>% filter(stop_sequence > 1 & type == "arrival") %>%
        mutate(segment_index = stop_sequence - 1)
    xFrom %>% select(route_id, trip_id, segment_index, time) %>%
        inner_join(
            xTo %>% select(segment_index, time), 
            by = "segment_index", 
            suffix = c("_from", "_to")
        ) %>% 
        mutate(
            time_from = time_from,
            time_to = time_to,
            travel_time = as.numeric(time_to - time_from)
        )
}
# x <- tu %>% filter(trip_id == tu$trip_id[1])
# x <- x %>% filter(vehicle_id == x$vehicle_id %>% table %>% sort %>% tail(1) %>% names)

tu <- tripupdates
rm(tripupdates)
if (file.exists("tu0.rda")) {
    load("tu0.rda")
} else {
    tu0 <- tu %>% distinct() %>%
        mutate(date = as.Date(format(timestamp, '%Y-%m-%d'))) %>%
        group_by(date, trip_id, vehicle_id) %>%
        do((.) %>% trip_times_to_travel_times())
    save(tu0, file = "tu0.rda")
}

library(RSQLite)
library(dbplyr)
rids <- unique(tu0$route_id)
con <- dbConnect(SQLite(), "fulldata.db")
routei <- con %>% tbl("routes") %>%
    filter(route_id %in% rids) %>%
    select(route_id, route_short_name, route_long_name) %>%
    collect()
dbDisconnect(con)
rnums <- structure(routei$route_short_name, .Names = routei$route_id)

tu1 <- tu0 %>% mutate(
    short_route_id = gsub('-.*', '', route_id),
    route_number = rnums[route_id]
)

    #     filter(between(
    #         timestamp, 
    #         median(timestamp) - 60 * 60, 
    #         median(timestamp) + 60 * 60
    #     ))
    # ) %>%
    # ungroup() %>%
    # group_by(vehicle_id, route_id, trip_id, timestamp) %>%
    # arrange(timestamp) %>%
    # summarize(
    #     stop_sequence = first(stop_sequence),
    #     type = first(type),
    #     time = first(time)
    # )

std <- function(x) {
    if (length(x) == 1) return(0)
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

rids <- names(sort(table(tu1$short_route_id), TRUE))
RID <- rids[21]; routei %>% filter(grepl(RID, route_id))

tu1r <- tu1 %>% 
    filter(
        short_route_id == RID &
        between(travel_time, 20, 60*30)
    ) %>%
    filter(
        segment_index < max((.)$segment_index)
    ) %>%
    group_by(segment_index) %>%
    do((.) %>% mutate(
        time = as.POSIXct(paste(Sys.Date(), format(time_from, '%H:%M:%S'))),
        dow = as.numeric(format(time_from, '%u')),
        week = as.factor(ifelse(as.numeric(format(time_from, '%d')) < 8, 1, 2)),
        weekend = as.numeric(dow > 5),
        travel_time.std = std(travel_time),
        tt_center = travel_time - mean(travel_time, na.rm = TRUE)
    )) %>%
    ungroup() %>% group_by(date, segment_index) %>%
    do((.) %>% arrange(time_from) %>% mutate(    
        delta = c(NA, diff(as.integer(time))),
        tt_diff = c(NA, diff(travel_time))
    )) %>%
    ungroup() %>%
    arrange(segment_index, time_from)
 
ggplot(
    tu1r, 
    aes(time, tt_center, colour = as.factor(week))
) + geom_point() + facet_grid(segment_index~weekend, scales = "free_y") +
    scale_x_datetime(label = function(x) format(x, '%H'))

ggplot(
    tu1r %>% filter(segment_index == 2), 
    aes(time, travel_time/60, colour = as.factor(dow))
) + geom_path(aes(group = date)) + facet_grid(week~weekend) +
    scale_x_datetime(label = function(x) format(x, '%H'))

ggplot(
    tu1r %>% filter(segment_index == 2), 
    aes(time, travel_time.std/60, colour = as.factor(dow))
) + geom_point() + facet_grid(week~weekend) +
    scale_x_datetime(label = function(x) format(x, '%H'))

ggplot(
    tu1r %>% filter(segment_index == 2),
    aes(time, tt_diff, colour = as.factor(dow))
) + geom_path(aes(group = dow)) + facet_grid(week ~ weekend) + 
    scale_x_datetime(label = function(x) format(x, '%H'))    

ggplot(
    tu1r %>% filter(segment_index == 2),
    aes(delta, tt_diff, colour = as.factor(dow))
) + geom_point() + facet_grid(week ~ weekend)

ggplot(
    tu1r,
    aes(time, abs(tt_diff / delta), colour = as.factor(segment_index))
) + geom_point() + facet_grid(week ~ weekend) +
    scale_x_datetime(label = function(x) format(x, '%H')) +
    geom_smooth(se=F)

ggplot(tu1r, aes(tt_diff / delta)) +
    geom_histogram() + facet_grid(week ~ weekend)


## generate a model to estimate the network system noise parameter, Q
library(rjags)
library(tidybayes)

tu1r1 <- tu1r %>% filter(week == 1) %>% arrange(segment_index, time_from)
jags.data <- 
    list(
        z = tu1r1$tt_center,
        delta = tu1r1$delta,
        s = tu1r1$segment_index,
        weekend = tu1r1$weekend,
        N = nrow(tu1r1),
        M = length(unique(tu1r1$segment_index)),
        J = length(unique(paste(tu1r1$segment_index, tu1r1$date))),
        Nj = tapply(1:nrow(tu1r1), paste(tu1r1$segment_index, tu1r1$date), min) %>% as.integer,
        Ni = tapply(1:nrow(tu1r1), paste(tu1r1$segment_index, tu1r1$date), length) %>% as.integer
    )

jags.fit <- 
    jags.model(
        'scripts/model1.jags', 
        jags.data,
        n.adapt = 1000,
        n.chains = 2
    )

samps <- 
    coda.samples(
        jags.fit, 
        c('x', 'r_sig', 'q', 'qsig', 'q_sig'), 
        n.iter = 10000, 
        thin = 10
    )

samps %>% spread_draws(q) %>%
    ggplot(aes(.iteration, group = .chain, colour = .chain %>% as.factor)) +
    geom_path(aes(y = q))

samps %>% spread_draws(qsig) %>%
    ggplot(aes(.iteration, group = .chain, colour = .chain %>% as.factor)) +
    geom_path(aes(y = qsig))

egg::ggarrange(
    samps %>% spread_draws(q_sig[segment_index]) %>%
        ggplot(aes(.iteration, q_sig, group = .chain)) + 
        facet_grid(segment_index~., scales = "free_y") + 
        geom_path(),
    samps %>% spread_draws(r_sig[segment_index]) %>%
        ggplot(aes(.iteration, r_sig, group = .chain)) + 
        facet_grid(segment_index~., scales = "free_y") + 
        geom_path()
)

p <- samps %>%
    spread_draws(r_sig[segment_index], q_sig[segment_index]) %>%
    # median_qi() %>%
    ggplot(aes(y = segment_index))

p + stat_pointintervalh(aes(x = r_sig))
p + stat_pointintervalh(aes(x = q_sig))

st <- samps %>% spread_draws(x[index])

tu1r1 %>% 
    bind_cols(
        st %>% filter(.iteration == sample(900:999, 1) & .chain == 1)
        # st %>% median_qi()
    ) %>%
    ggplot(aes(time, tt_center, colour = as.factor(dow))) + 
    geom_point() + 
    geom_path(aes(time, x)) +
    facet_grid(segment_index~weekend, scales = "free_y") +
    scale_x_datetime(label = function(x) format(x, '%H'))


## segment lengths
con <- dbConnect(SQLite(), "fulldata.db")
tids <- unique(tu1r1$trip_id)
tsegs <- con %>% tbl("trips") %>% 
    filter(trip_id %in% tids) %>%
    left_join(con %>% tbl("shape_segments")) %>%
    left_join(con %>% tbl("road_segments") %>% select(road_segment_id, length)) %>%
    select(trip_id, road_segment_id, shape_road_sequence, length) %>%
    collect()
dbDisconnect(con)

tu1r1 %>% left_join(tsegs %>% select(trip_id, shape_road_sequence, length), 
    by = c("trip_id", "segment_index" = "shape_road_sequence")) %>%
    select(segment_index, length, travel_time) %>% distinct %>%
    left_join(samps %>% spread_draws(q_sig[segment_index]) %>% median_qi()) %>%
    ggplot(aes(travel_time, q_sig)) + geom_point()

tu1r1 %>% left_join(tsegs %>% select(trip_id, shape_road_sequence, length), 
    by = c("trip_id", "segment_index" = "shape_road_sequence")) %>%
    select(segment_index, length, travel_time) %>% distinct %>%
    left_join(samps %>% spread_draws(r_sig[segment_index]) %>% median_qi()) %>%
    ggplot(aes(travel_time, r_sig)) + geom_point()

qr.ests <- samps %>% spread_draws(q_sig[seg], r_sig[seg]) %>% median_qi()

## apply KF to second week with values of R and Q and calculate RMSE
tu1r2 <- tu1r %>% filter(week == 2) %>% arrange(segment_index, time_from)
test.data <- tu1r2 %>% 
    mutate(
        z = tt_center,
        t = time,
        s = segment_index,
        ind = as.integer(as.factor(paste(segment_index, date))),
        i = 1:n()
    ) %>%
    select(z, t, delta, s, weekend, ind, i)

KF <- function(d, x0, p0, q, r) {
    # ggplot(d, aes(t, z)) + geom_point() + 
    #     geom_pointrange(aes(min(d$t), x0, ymin = x0 - p0, ymax = x0 + p0), 
    #         data=NULL, col="red")
    xhat <- numeric(nrow(d))
    xhat.pred <- numeric(nrow(d))
    p <- numeric(nrow(d))
    p.pred <- numeric(nrow(d))

    Q <- q^2
    R <- r^2

    for (i in seq_along(1:nrow(d))) {
        # predict
        if (i == 1) {
            xhat.pred[i] <- x0
            p.pred[i] <- p0
        } else {
            xhat.pred[i] <- xhat[i-1]
            p.pred[i] <- p[i-1] + Q * d$delta[i]
        }

        # update
        y <- d$z[i] - xhat.pred[i]
        S <- R + p.pred[i]
        K <- p.pred[i] / S
        xhat[i] <- xhat.pred[i] + K * y
        p[i] <- (1 - K) * p.pred[i] * (1 - K) + K * R * K
    }
    d <- d %>% mutate(
        xhat = xhat,
        xhat.pred = xhat.pred,
        p = p,
        p.pred = p.pred
        )

    # ggplot(d, aes(t)) +
    #     geom_point(aes(y = z)) +
    #     geom_pointrange(aes(y = xhat, ymin = xhat - p, ymax = xhat + p), colour = "blue") +
    #     geom_pointrange(aes(y = xhat.pred, ymin = xhat.pred - p.pred, ymax = xhat.pred + p.pred), 
    #         colour = "red")
    d
}

test.results <- test.data %>% group_by(ind) %>%
    do({
        d <- (.)
        x0 <- st %>% filter(index == d$i[1]) %>% pluck("x") %>% mean
        p0 <- st %>% filter(index == d$i[1]) %>% pluck("x") %>% sd
        q <- qr.ests %>% filter(seg == d$s[1]) %>% pluck("q_sig")
        r <- qr.ests %>% filter(seg == d$s[1]) %>% pluck("r_sig")
        KF(d, x0, p0, q, r)
    })

test.results %>% 
    ggplot(aes(t, z, group = ind)) + 
    geom_point() + 
    geom_path(aes(t, xhat), col = "blue") +
    geom_path(aes(t, xhat.pred), col = "red") +
    facet_grid(s~weekend, scales = "free_y") +
    scale_x_datetime(label = function(x) format(x, '%H'))

# rmse
with(test.results, z - xhat)^2 %>% mean
with(test.results, z - xhat.pred)^2 %>% mean


########### Model 2 (2nd order)

jags.fit <- 
    jags.model(
        'scripts/model2.jags', 
        jags.data,
        n.adapt = 1000,
        n.chains = 2
    )

samps <- 
    coda.samples(
        jags.fit, 
        c('x', 'xdot', 'r_sig', 'q_sig', 'qw_sig'), 
        n.iter = 10000, 
        thin = 10
    )

# samps %>% spread_draws(q) %>%
#     ggplot(aes(.iteration, group = .chain, colour = .chain %>% as.factor)) +
#     geom_path(aes(y = q))

# samps %>% spread_draws(qsig) %>%
#     ggplot(aes(.iteration, group = .chain, colour = .chain %>% as.factor)) +
#     geom_path(aes(y = qsig))

egg::ggarrange(
    nrow = 1,
    samps %>% spread_draws(q_sig[segment_index]) %>%
        ggplot(aes(.iteration, q_sig, group = .chain)) + 
        facet_grid(segment_index~., scales = "free_y") + 
        geom_path(),
    samps %>% spread_draws(qw_sig[segment_index]) %>%
        ggplot(aes(.iteration, qw_sig, group = .chain)) + 
        facet_grid(segment_index~., scales = "free_y") + 
        geom_path(),
    samps %>% spread_draws(r_sig[segment_index]) %>%
        ggplot(aes(.iteration, r_sig, group = .chain)) + 
        facet_grid(segment_index~., scales = "free_y") + 
        geom_path()
)

p <- samps %>%
    spread_draws(r_sig[segment_index], q_sig[segment_index]) %>%
    # median_qi() %>%
    ggplot(aes(y = segment_index))

p + stat_pointintervalh(aes(x = r_sig))
p + stat_pointintervalh(aes(x = q_sig))

st <- samps %>% spread_draws(x[index], xdot[index])

pxx <- tu1r1 %>% 
    bind_cols(
        st %>% filter(.iteration == sample(900:999, 1) & .chain == 1)
    ) %>% 
    ggplot(aes(time, tt_center, colour = as.factor(dow))) + 
    facet_grid(segment_index~weekend, scales = "free_y") +
    scale_x_datetime(label = function(x) format(x, '%H'))

egg::ggarrange(
    nrow = 1,
    pxx + geom_point() + geom_path(aes(time, x)),
    pxx + geom_path(aes(time, xdot))
)


## segment lengths
con <- dbConnect(SQLite(), "fulldata.db")
tids <- unique(tu1r1$trip_id)
tsegs <- con %>% tbl("trips") %>% 
    filter(trip_id %in% tids) %>%
    left_join(con %>% tbl("shape_segments")) %>%
    left_join(con %>% tbl("road_segments") %>% select(road_segment_id, length)) %>%
    select(trip_id, road_segment_id, shape_road_sequence, length) %>%
    collect()
dbDisconnect(con)

tu1r1 %>% left_join(tsegs %>% select(trip_id, shape_road_sequence, length), 
    by = c("trip_id", "segment_index" = "shape_road_sequence")) %>%
    select(segment_index, length, travel_time) %>% distinct %>%
    left_join(samps %>% spread_draws(q_sig[segment_index]) %>% median_qi()) %>%
    ggplot(aes(travel_time, q_sig)) + geom_point()

tu1r1 %>% left_join(tsegs %>% select(trip_id, shape_road_sequence, length), 
    by = c("trip_id", "segment_index" = "shape_road_sequence")) %>%
    select(segment_index, length, travel_time) %>% distinct %>%
    left_join(samps %>% spread_draws(r_sig[segment_index]) %>% median_qi()) %>%
    ggplot(aes(travel_time, r_sig)) + geom_point()

qr.ests <- samps %>% spread_draws(q_sig[seg], r_sig[seg]) %>% median_qi()

## apply KF to second week with values of R and Q and calculate RMSE
tu1r2 <- tu1r %>% filter(week == 2) %>% arrange(segment_index, time_from)
test.data2 <- tu1r2 %>% 
    mutate(
        z = tt_center,
        t = time,
        s = segment_index,
        ind = as.integer(as.factor(paste(segment_index, date))),
        i = 1:n()
    ) %>%
    select(z, t, delta, s, weekend, ind, i)

KF2 <- function(d, x0, p0, q, r) {
    # ggplot(d, aes(t, z)) + geom_point() + 
    #     geom_pointrange(aes(min(d$t), x0, ymin = x0 - p0, ymax = x0 + p0), 
    #         data=NULL, col="red")
    xhat <- matrix(nrow = nrow(d), ncol = 2)
    xhat.pred <- matrix(nrow = nrow(d), ncol = 2)
    p <- matrix(nrow = nrow(d), ncol = 4)
    p.pred <- matrix(nrow = nrow(d), ncol = 4)

    Q <- diag(c(0, q^2))
    R <- r^2
    F <- function(delta) matrix(c(1, 0, ifelse(is.na(delta), 0, delta), 1), nrow = 2)
    H <- rbind(c(1, 0))
    I <- diag(2)

    for (i in seq_along(1:nrow(d))) {
        # predict
        Fi <- F(d$delta[i])
        if (i == 1) {
            xhat.pred[i,] <- x0
            p.pred[i,] <- p0
            X <- t(xhat.pred[i, , drop = FALSE])
            P <- matrix(p.pred[i, ], nrow = 2)
        } else {
            X <- Fi %*% X
            P <- Fi %*% P %*% t(Fi) + Q
            xhat.pred[i, ] <- c(X)
            p.pred[i, ] <- c(P)
        }

        # update
        y <- d$z[i] - H %*% X
        S <- R + H %*% P %*% t(H)
        K <- P %*% t(H) %*% solve(S)
        X <- X + K %*% y
        P <- (I - K %*% H) * P * t(I - K %*% H) + K %*% R %*% t(K)
        xhat[i, ] <- c(X)
        p[i, ] <- c(P)
    }
    d <- d %>% 
        mutate(
            x_hat = xhat[,1],
            xd_hat = xhat[,2],
            x_hat.pred = xhat.pred[,1],
            xd_hat.pred = xhat.pred[,2],
            p11 = p[,1],
            p21 = p[,2],
            p12 = p[,3],
            p22 = p[,4],
            p.pred11 = p.pred[,1],
            p.pred21 = p.pred[,2],
            p.pred12 = p.pred[,3],
            p.pred22 = p.pred[,4]
        )

    # ggplot(d, aes(t)) +
    #     geom_point(aes(y = z)) +
    #     geom_pointrange(aes(y = xhat, ymin = xhat - p, ymax = xhat + p), colour = "blue") +
    #     geom_pointrange(aes(y = xhat.pred, ymin = xhat.pred - p.pred, ymax = xhat.pred + p.pred), 
    #         colour = "red")
    d
}

test.results2 <- test.data2 %>% group_by(ind) %>%
    do({
        d <- (.)
        x0 <- st %>% filter(index == d$i[1]) %>% ungroup %>% select(x, xdot) %>% colMeans
        p0 <- st %>% filter(index == d$i[1]) %>% ungroup %>% select(x, xdot) %>% cov %>% c
        q <- qr.ests %>% filter(seg == d$s[1]) %>% pluck("q_sig")
        r <- qr.ests %>% filter(seg == d$s[1]) %>% pluck("r_sig")
        KF2(d, x0, p0, q, r)
    })

ptt <- test.results2 %>% 
    ggplot(aes(t, z, group = ind)) + 
    facet_grid(s~weekend, scales = "free_y") +
    scale_x_datetime(label = function(x) format(x, '%H'))

egg::ggarrange(
    nrow = 1,
    ptt + geom_point() + 
        geom_path(aes(t, x_hat), col = "blue") +
        geom_path(aes(t, x_hat.pred), col = "red"),
    ptt + geom_path(aes(t, xd_hat), col = "blue") +
        geom_path(aes(t, xd_hat.pred), col = "red")
)

# rmse
with(test.results2, z - x_hat)^2 %>% mean
with(test.results2, z - x_hat.pred)^2 %>% mean


pt1 <- test.results %>% 
    ggplot(aes(t, group = ind)) + 
    facet_grid(s~weekend, scales = "free_y") +
    scale_x_datetime(label = function(x) format(x, '%H'))
pt2 <- test.results2 %>%
    ggplot(aes(t, group = ind)) +
    facet_grid(s~weekend, scales = "free_y") +
    scale_x_datetime(label = function(x) format(x, '%H'))

pt1 + geom_linerange(aes(ymin = 0, ymax = xhat - xhat.pred))
pt2 + geom_linerange(aes(ymin = 0, ymax = x_hat - x_hat.pred))

pt1 + geom_linerange(aes(ymin = 0, ymax = z - xhat.pred))
pt2 + geom_linerange(aes(ymin = 0, ymax = z - x_hat.pred))

pt1 + geom_linerange(aes(ymin = 0, ymax = z - xhat))
pt2 + geom_linerange(aes(ymin = 0, ymax = z - x_hat))

test.results %>% 
    ungroup() %>%
    mutate(
        err_obs_pred = abs(z - xhat.pred),
        err_obs_est = abs(z - xhat),
        err_est_pred = abs(xhat - xhat.pred)
    ) %>%
    gather(key = "type", value = "error", err_obs_pred, err_obs_est, err_est_pred) %>%
    select(error, type) %>%
    mutate(test = "test1") %>%
    bind_rows(
        test.results2 %>% 
            ungroup() %>%
            mutate(
                err_obs_pred = abs(z - x_hat.pred),
                err_obs_est = abs(z - x_hat),
                err_est_pred = abs(x_hat - x_hat.pred)
            ) %>%
            gather(key = "type", value = "error", err_obs_pred, err_obs_est, err_est_pred) %>%
            select(error, type) %>%
            mutate(test = "test2")
    ) %>%
    ggplot(aes(test, error, fill = test)) +
    facet_grid(~type) +
    geom_violin() 



## num -> [30min intervals]
times <- pretty(tu1r$time, 20 * 2)
smry <- tu1r %>%
    mutate(
        period = cut((.)$time, breaks = times)
    ) %>%
    group_by(date, segment_index, period) %>%
    summarize(
        tt = mean(travel_time, na.rm = TRUE),
        sd.tt = sd(travel_time, na.rm = TRUE),
        n = n(),
        dow = first(dow),
        week = first(week),
        weekend = first(weekend)
    ) %>%
    mutate(
        time = as.POSIXct(period)
    )

ggplot(smry %>% filter(segment_index == 1), 
    aes(time, tt/60, colour = as.factor(dow))) + 
    geom_point() +
    facet_grid(week~weekend) +
    scale_x_datetime(label = function(x) format(x, '%H')) 


## drop it into a matrix




## now a combination of [trip_id, shape_id, segment_index] can get [segment_id]
con <- dbConnect(SQLite(), "fulldata.db")
tids <- unique(tu1$trip_id)
tsegs <- con %>% tbl("trips") %>% 
    filter(trip_id %in% tids) %>%
    left_join(con %>% tbl("shape_segments")) %>%
    left_join(con %>% tbl("road_segments") %>% select(road_segment_id, length)) %>%
    select(trip_id, road_segment_id, shape_road_sequence, length) %>%
    collect()
segids <- unique(tsegs$road_segment_id)
segs <- con %>% tbl("road_segments") %>%
    filter(road_segment_id %in% segids) %>%
    left_join(con %>% tbl("intersections"), by = c("int_from" = "intersection_id")) %>%
    left_join(con %>% tbl("intersections"), by = c("int_to" = "intersection_id"), suffix = c("", "_to")) %>%
    collect()
dbDisconnect(con)

## short trip ...
tu1 <- tu1 %>% mutate(short_trip_id = gsub('-.*', '', trip_id))
tsegs <- tsegs %>% mutate(short_trip_id = gsub('-.*', '', trip_id))
sx <- tu1 %>% inner_join(tsegs, by = c("short_trip_id", "segment_index" = "shape_road_sequence")) %>%
    filter(travel_time > 20 & travel_time < 60*30) %>%
    arrange(road_segment_id, time_from) %>%
    mutate(speed = length / travel_time) %>%
    filter(speed < 30)

NX1segs <- sx %>% ungroup() %>% 
    filter(grepl(RID, route_id)) %>%
    select(road_segment_id, segment_index) %>% distinct() %>%
    arrange(segment_index) %>%
    filter(segment_index < max(segment_index))

# seg20 <- sx$road_segment_id %>% table %>% sort(TRUE) %>% names
sxnex <- sx %>% filter(road_segment_id %in% NX1segs$road_segment_id) %>%
    mutate(
        segment = factor(road_segment_id, levels = NX1segs$road_segment_id),
        time = as.POSIXct(paste(Sys.Date(), format(time_from, '%H:%M:%S'))),
        dow = as.numeric(format(time_from, '%u')),
        week = as.factor(ifelse(as.numeric(format(time_from, '%d')) < 8, 1, 2)),
        weekend = as.numeric(dow > 5)
    )

times2 <- pretty(sxnex$time, 20 * 2)
smry2 <- sxnex %>% ungroup () %>%
    mutate(
        period = cut((.)$time, breaks = times2)
    ) %>%
    group_by(date, segment, period) %>%
    summarize(
        tt = mean(travel_time, na.rm = TRUE),
        sd.tt = sd(travel_time, na.rm = TRUE),
        n = n(),
        speed = mean(speed, na.rm = TRUE)
    ) %>%
    ungroup () %>%
    mutate(
        time = as.POSIXct(period),
        dow = as.numeric(format(date, '%u')),
        week = as.factor(ifelse(as.numeric(format(date, '%d')) < 8, 1, 2)),
        weekend = as.numeric(dow > 5)
    )


ggplot(smry2, aes(time, speed / 1000 * 60 * 60)) + 
# ggplot(smry2, aes(time, tt / 60)) + 
    geom_point(aes(colour = as.factor(week))) +
    # geom_path(aes(colour = date)) +
    facet_grid(segment~weekend, scales = "free_y") +
    xlab("Time") + 
    ylab("Travel time (min)") +
    scale_x_datetime(label = function(x) format(x, '%H:%M'))
    # geom_smooth(aes(colour = route_number), se = FALSE) + 
    # scale_colour_discrete(guide = FALSE) + 
    # scale_y_log10()
    # geom_smooth(method = lm, formula = y ~ splines::bs(x, 4))


### I want the model to work on both mu and sigma
# i.e., mu(t) and sigma(t)

alldata <- smry2 %>% 
    select(date, time, segment, tt, speed, dow, week, weekend) %>%
    mutate(
        date = factor(date),
        dow = factor(dow),
        dayofweek = factor(format(as.Date(date), '%a'), 
            levels = c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun')),
        week = factor(week),
        saturday = as.numeric(dow == 6),
        sunday = as.numeric(dow == 7)
    ) #%>% filter(segment == '4691')
train <- alldata %>% filter(week == 1)
test <- alldata %>% filter(week == 2)

rmse <- function(fit, data) 
    mean((predict(fit, newdata = data) - data$tt)^2, na.rm = TRUE)


fit0 <- lm(tt ~ 1, data = train)
rmse(fit0, test)

fit1 <- lm(tt ~ segment * dow - 1, data = train)
rmse(fit1, test)

fit2 <- lm(tt ~ segment + segment*saturday + segment*sunday - 1, data = train)
rmse(fit2, test)


## spline models
library(splines)
fitSpline <- function(i, j)
    lm(speed ~ bs(time, intercept = TRUE, degree = i):segment + 
        bs(time, intercept = TRUE, degree = j):weekend:segment - 1, data = train)

RMSE <- expand.grid(i = 1:20, j = 1:15)
RMSE$rmse <- apply(RMSE, 1, function(x) {
    rmse(fitSpline(x["i"], x["j"]), test)
})
ggplot(RMSE, aes(i, j, size = log(rmse))) + geom_point()

RMSE[which.min(RMSE$rmse), ]
sfit <- fitSpline(5, 2)
summary(sfit)
test %>% mutate(speed.pred = predict(sfit, newdata = test)) %>%
    ggplot(aes(time, speed.pred*3.6)) + 
    geom_point(aes(y = speed*3.6), data = train, col = 'gray', size = 0.5) +
    geom_point(aes(y = speed*3.6), data = test, col = 'orangered', size = 0.5) +
    geom_path(aes(group = date), colour = 'steelblue') +
    facet_grid(segment~dayofweek, scales = "free_y") +
    xlab("Time") + 
    ylab("Speed (m/s)") + #ylim(0, 100) +
    scale_x_datetime(label = function(x) format(x, '%H:%M'))
    # ylab("Travel time (min)") +


## lag model
alldata <- alldata %>% mutate(time_index = time %>% as.factor %>% as.integer)

alldata.lag <- alldata %>%
    filter(time_index > 1) %>%
    mutate(time_index_lag = time_index - 1) %>%
    left_join(
        alldata %>% select(date, segment, time_index, tt, speed),
        by = c('date', 'segment', 'time_index_lag' = 'time_index'),
        suffix = c("", ".lag")
    )


train.lag <- alldata.lag %>% filter(week == 1)
test.lag <- alldata.lag %>% filter(week == 2)

ggplot(train.lag, aes(time, tt - tt.lag)) +
    geom_point() +
    facet_grid(segment~dayofweek, scales = "free_y") +
    scale_x_datetime(label = function(x) format(x, '%H:%M'))


fitLagSpline <- function(i, j)
    lm(I(tt - tt.lag) ~ bs(time, degree = i, intercept = TRUE):segment +
        bs(time, degree = j, intercept = TRUE):weekend:segment - 1, data = train.lag)

RMSE <- expand.grid(i = 1:10, j = 1:5)
RMSE$rmse <- apply(RMSE, 1, function(x) {
    rmse(fitLagSpline(x["i"], x["j"]), test.lag)
})
ggplot(RMSE, aes(i, j, size = rmse)) + geom_point()

RMSE[which.min(RMSE$rmse), ]
lfit <- fitLagSpline(8, 3)
summary(lfit)

test.lag %>% mutate(tt.pred = predict(lfit, newdata = test.lag)) %>%
    ggplot(aes(time, (tt.pred)/60)) + 
    geom_point(aes(y = (tt - tt.lag)/60), data = train.lag, col = 'gray', size = 0.5) +
    geom_path(aes(y = (tt.pred)/60, group = date), 
        data = train.lag %>% mutate(tt.pred = predict(lfit, train.lag)), col = 'gray') +
    geom_point(aes(y = (tt - tt.lag)/60), col = 'orangered', size = 0.5) +
    # geom_path(col = 'orangered') +
    facet_grid(segment~weekend, scales = "free_y") +
    xlab("Time") + 
    ylab("Travel Time (min)") + 
    scale_x_datetime(label = function(x) format(x, '%H:%M'))
    # ggplot(aes(time, speed.pred*3.6)) + 
    # geom_point(aes(y = speed*3.6), data = train, col = 'gray', size = 0.5) +
    # geom_point(aes(y = speed*3.6), data = test, col = 'orangered', size = 0.5) +
    # ylab("Travel time (min)") +


###### NN model ...
scaleX <- function(x) (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE)) * 2 - 1
nn.data <- alldata %>%
    inner_join(NX1segs %>% mutate(segment = as.factor(road_segment_id))) %>%
    group_by(date, segment_index) %>%
    do((.) %>% arrange(time) %>% mutate(speed_lag = c(NA, speed[-n()]))) %>%
    ungroup() %>% group_by(date, time) %>%
    do(
        (.) %>% arrange(segment_index) %>%
        mutate(
            speed_prev_seg = c(NA, speed_lag[-n()]),
            speed_next_seg = c(speed_lag[-1], NA)
        )
    ) %>%
    ungroup() %>%
    mutate(
        time = scaleX(as.integer(time)),
        speed = scaleX(speed),
        speed_lag = scaleX(speed_lag),
        speed_prev_seg = scaleX(speed_prev_seg),
        speed_next_seg = scaleX(speed_next_seg)
    )

nn.data %>% 
    ggplot(aes(time, speed, colour = week)) +
    geom_point() +
    facet_wrap(~segment_index)

nn.data %>% 
    ggplot(aes(speed, speed_next_seg, colour = dow)) + 
    geom_point() +
    facet_wrap(~segment_index, scales="free")

nn.train <- nn.data %>% filter(week == 1)
nn.test <- nn.data %>% filter(week == 2)

nn2 <- nn.train %>% filter(!is.na(speed_lag) & !is.na(speed_prev_seg) & !is.na(speed_next_seg))
library(neuralnet)
fit.nn <- neuralnet(speed ~ time + weekend + #as.factor(segment_index) + dayofweek +
    speed_lag + speed_prev_seg + speed_next_seg,
    data = nn2, hidden = c(3), linear.output = TRUE)
plot(fit.nn)

nn.test %>% 
    bind_cols(
        speed_pred = compute(
            fit.nn, 
            nn.test %>% select(time, weekend, speed_lag, speed_prev_seg, speed_next_seg)
        )$net.result
    ) %>%
    ggplot(aes(time)) + 
    geom_point(aes(speed)) +
    geom_point(aes(y = speed_pred), colour = "red")
