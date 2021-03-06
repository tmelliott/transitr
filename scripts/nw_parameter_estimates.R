library(tidyverse)
library(ggraph)
library(tidygraph)
library(Matrix)
library(rjags)
library(RSQLite)
library(dbplyr)

source("scripts/common2.R")

sim <- "sim000"

### --- work on a SUBSET (route NX1)

## fetch the data
segdata <- get_segment_data("982")
nwobs <- 
    read_csv(
        file.path("simulations", sim, "segment_observations.csv"),
        col_names = c("segment_id", "timestamp", "travel_time", "measurement_error"),
        col_types = "iiin"
    ) %>% 
    filter(timestamp > 0 & segment_id %in% unique(segdata$road_segment_id)) %>%
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

sids <- unique(segdata$road_segment_id)

con <- dbConnect(SQLite(), "fulldata.db")
seglens <- con %>% tbl("road_segments") %>% 
    filter(road_segment_id %in% sids) %>%
    select(road_segment_id, length, int_from, int_to) %>% collect %>%
    mutate(road_segment_id,
           from = int_from, to = int_to)
dbDisconnect(con)

nwobs <- nwobs %>% 
    left_join(seglens, by = c("segment_id" = "road_segment_id")) %>%
    mutate(length = length.x)
g <- nwobs %>%
    group_by(segment_id) %>% 
        summarize(
            speed = mean(speed) * 60 * 60 / 1000,
            from = first(from), to = first(to)
        ) %>%
    as_tbl_graph()

ggraph(g, layout = 'kk') + 
    geom_edge_fan(
        aes(color = speed),
        arrow = arrow(length = unit(1, 'mm')), 
        start_cap = circle(1, 'mm'),
        end_cap = circle(1, 'mm')) +
    scale_edge_colour_gradientn(colours = viridis::viridis(256, option = "A"), limits = c(0, 100)) +
    geom_node_point(shape = 19, size = 0.5)



## look at the data
seg <- nwobs$segment_id %>% table %>% sort %>% names %>% tail(10)

ps <- ggplot(nwobs) #+ #%>% filter(segment_id %in% seg)) +

# ggplot(nwobs, aes(timestamp, speedkm)) + 
#     geom_point(alpha = 0.05) +
#     geom_quantile(quantiles = 0.5)

ps + geom_histogram(aes(travel_time)) + scale_x_continuous(trans = "log")
ps + geom_histogram(aes(speedkm)) + xlim(0, 100)

ps + geom_path(aes(timestamp, travel_time)) + 
    facet_wrap(~segment_id, scales="free_y") +
    scale_x_datetime(labels = function(x) 
        paste0(as.numeric(format(x, "%I")), format(x, "%P")))
ps + geom_point(aes(timestamp, speedkm)) + ylim(0, 100)

## which time period does each obs belong?
times <- pretty(nwobs$timestamp, 20 * 2)
nwobsf <- nwobs %>%
    # filter(segment_id %in% seg) %>%
    mutate(period = cut((.)$timestamp, breaks = times)) %>%
    group_by(segment_id, period) %>%
    do(
        (.) %>% mutate(wt = 1 / sum(1 / measurement_error^2) / measurement_error^2)
    ) %>%
    summarize(
        sumwt = sum(wt),
        time = sum(travel_time * wt, na.rm = TRUE),
        error = sum(wt * (travel_time - time)^2),
        error = sd(travel_time, na.rm = TRUE),
        # log_time = mean(log(travel_time), na.rm = TRUE),
        # log_error = sd(log(travel_time), na.rm = TRUE),
        n = n(),
        length = first(length)
    ) %>%
    ungroup () %>%
    mutate(
        speed = length / exp(time),
        speedkm = speed / 1000 * 60 * 60
    ) %>%
    arrange(segment_id, period) %>%
    mutate(
        period_int = as.integer(period),
        period = as.POSIXct(period)
    )

Ymat <- nwobsf %>% 
    select(segment_id, period_int, time) %>% 
    spread(key = period_int, value = time) %>%
    select(-1) %>% as.matrix
Yerr <- nwobsf %>% 
    select(segment_id, period_int, error) %>%
    spread(key = period_int, value = error) %>%
    select(-1) %>% as.matrix

ggplot(nwobsf) +
    geom_ribbon(aes(period, ymin = pmax(0, time - error), ymax = time + error),
        fill = "gray") +
    geom_path(aes(period, time)) +
    facet_wrap(~segment_id, scales = "free_y") +
    scale_x_datetime(labels = function(x) 
        paste0(as.numeric(format(x, "%I")), format(x, "%P")))

mean(apply(Ymat, 1, function(x) sd(diff(x), na.rm = TRUE)), na.rm = TRUE) ## Sigma_mu
sd(apply(Ymat, 1, function(x) sd(diff(x), na.rm = TRUE)), na.rm = TRUE)   ## Sigma_sig

nwdata <- 
    list(
        M = nrow(Ymat),
        T = ncol(Ymat),
        Y = Ymat %>% log
    )

jags.fit <- 
    jags.model(
        'scripts/nw_model_indep.jags', 
        nwdata,
        n.adapt = 10000,
        n.chains = 1
    )

samps <- 
    coda.samples(
        jags.fit, 
        c('Q_mu', 'Q_sig', 'R_mu', 'R_sig'), 
        n.iter = 10000, 
        thin = 10
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
    as.tibble() %>%
    arrange(time, segment)
i <- sample(990, 10)
X1df <- xdf %>% 
    mutate(
        obs = as.numeric(Ymat),
        estimate0 = as.array(Xs)[i[10],] %>% matrix(nrow = nwdata$M) %>% as.numeric,
        estimate1 = as.array(Xs)[i[1],] %>% matrix(nrow = nwdata$M) %>% as.numeric,
        estimate2 = as.array(Xs)[i[2],] %>% matrix(nrow = nwdata$M) %>% as.numeric,
        estimate3 = as.array(Xs)[i[3],] %>% matrix(nrow = nwdata$M) %>% as.numeric,
        estimate4 = as.array(Xs)[i[4],] %>% matrix(nrow = nwdata$M) %>% as.numeric,
        estimate5 = as.array(Xs)[i[5],] %>% matrix(nrow = nwdata$M) %>% as.numeric,
        estimate6 = as.array(Xs)[i[6],] %>% matrix(nrow = nwdata$M) %>% as.numeric,
        estimate7 = as.array(Xs)[i[7],] %>% matrix(nrow = nwdata$M) %>% as.numeric,
        estimate8 = as.array(Xs)[i[8],] %>% matrix(nrow = nwdata$M) %>% as.numeric,
        estimate9 = as.array(Xs)[i[9],] %>% matrix(nrow = nwdata$M) %>% as.numeric
    )
ggplot(X1df) + 
    geom_point(aes(time, obs)) + 
    geom_path(aes(time, estimate0), colour = 'red') + 
    geom_path(aes(time, estimate1), colour = 'red') + 
    geom_path(aes(time, estimate2), colour = 'red') + 
    geom_path(aes(time, estimate3), colour = 'red') + 
    geom_path(aes(time, estimate4), colour = 'red') + 
    geom_path(aes(time, estimate5), colour = 'red') + 
    geom_path(aes(time, estimate6), colour = 'red') + 
    geom_path(aes(time, estimate7), colour = 'red') + 
    geom_path(aes(time, estimate8), colour = 'red') + 
    geom_path(aes(time, estimate9), colour = 'red') + 
    facet_wrap(~segment)

