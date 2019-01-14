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
seg <- nwobs$segment_id %>% table %>% sort %>% names %>% tail(10)

ps <- ggplot(nwobs %>% filter(segment_id %in% seg)) +
    facet_wrap(~segment_id)

# ggplot(nwobs, aes(timestamp, speedkm)) + 
#     geom_point(alpha = 0.05) +
#     geom_quantile(quantiles = 0.5)

ps + geom_histogram(aes(travel_time)) + scale_x_continuous(trans = "log")
ps + geom_histogram(aes(speedkm)) + xlim(0, 100)

ps + geom_path(aes(timestamp, travel_time)) #+ scale_y_log10()
ps + geom_point(aes(timestamp, speedkm)) + ylim(0, 100)


## which time period does each obs belong?
times <- pretty(nwobs$timestamp, 20 * 4)
nwobsf <- nwobs %>%
    filter(segment_id %in% seg) %>%
    mutate(period = cut((.)$timestamp, breaks = times)) %>%
    group_by(segment_id, period) %>%
    summarize(
        time = mean(travel_time, na.rm = TRUE),
        error = sd(travel_time, na.rm = TRUE),
        log_time = mean(log(travel_time), na.rm = TRUE),
        log_error = sd(log(travel_time), na.rm = TRUE),
        n = n(),
        length = first(length)
    ) %>%
    ungroup () %>%
    mutate(
        speed = length / exp(time),
        speedkm = speed / 1000 * 60 * 60
    ) %>%
    arrange(segment_id, period) %>%
    mutate(period_int = as.integer(period))

Ymat <- nwobsf %>% 
    select(segment_id, period_int, time) %>% 
    spread(key = period_int, value = time) %>%
    select(-1) %>% as.matrix
Yerr <- nwobsf %>% 
    select(segment_id, period_int, error) %>%
    spread(key = period_int, value = error) %>%
    select(-1) %>% as.matrix

ggplot(nwobsf) +
    geom_path(aes(period_int, time)) +
    facet_wrap(~segment_id)

mean(apply(Ymat, 1, function(x) sd(diff(x), na.rm = TRUE)), na.rm = TRUE) ## Sigma_mu
sd(apply(Ymat, 1, function(x) sd(diff(x), na.rm = TRUE)), na.rm = TRUE)   ## Sigma_sig

nwdata <- 
    list(
        M = nrow(Ymat),
        T = ncol(Ymat),
        Y = Ymat
    )

jags.fit <- 
    jags.model(
        'scripts/nw_model_indep.jags', 
        nwdata,
        n.adapt = 1000,
        n.chains = 1
    )

samps <- 
    coda.samples(
        jags.fit, 
        c('Q_mu', 'Q_sig', 'R_mu', 'R_sig'), 
        n.iter = 1000, 
        thin = 1
    )

summary(samps)
plot(samps)

# sds <- coda.samples(jags.fit, 'sig_r', n.iter = 10000, thin = 10)
# summary(sds)

# Sigmas <- coda.samples(jags.fit, 'Sigma', n.iter = 10000, thin = 10)
# summary(Sigmas)


Xs <- coda.samples(jags.fit, 'X', n.iter = 1000, thin = 1)

xdf <- 
    expand.grid(
        segment = 1:nwdata$M, 
        time = as.POSIXct(as.character(unique(nwobsf$period)), origin = "1970-01-01")
    ) %>% 
    as.tibble()
for (i in sample(1:10)) {
    X1 <- matrix(as.array(Xs)[i,], nrow = nwdata$M)
    X1df <- xdf %>% mutate(estimate = as.numeric(X1), obs = as.numeric(Ymat))

    g <- ggplot(X1df) + 
        geom_point(aes(time, obs)) + 
        geom_path(aes(time, estimate), colour = 'red') + 
        facet_wrap(~segment)
    dev.hold()
    print(g)
    dev.flush()
}

Xds <- coda.samples(jags.fit, 'Xdot', n.iter = 1000, thin = 1)

xdf <- 
    expand.grid(
        segment = 1:nwdata$M, 
        time = as.POSIXct(as.character(unique(nwobsf$period)), origin = "1970-01-01")
    ) %>% 
    as.tibble()
for (i in sample(1:10)) {
    Xd1 <- matrix(as.array(Xds)[i,], nrow = nwdata$M)
    Xd1df <- xdf %>% mutate(estimate = as.numeric(Xd1))

    g <- ggplot(Xd1df) + 
        geom_path(aes(time, estimate), colour = 'red') + 
        facet_wrap(~segment) + ylim(range(as.matrix(Xds)))
    dev.hold()
    print(g)
    dev.flush()
}
