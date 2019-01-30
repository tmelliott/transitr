library(tidyverse)

get_sim_files <- function(sim) {
    if (file.exists(file.path("simulations", sim, "modeleval.rds"))) {
        return(readRDS(file.path("simulations", sim, "modeleval.rds")))
    }
    x <- try({
        siminfo <- strsplit(sim, "_")[[1]][-1]
        if (grepl("e", siminfo[3])) siminfo[3] <- format(as.numeric(siminfo[3]), scientific = FALSE)
        siminfo <- as.numeric(gsub("-", ".", siminfo))
        do.call(bind_rows, 
            lapply(list.files(file.path("simulations", sim, "modeleval"), pattern="vehicle_.*\\.csv", full.names = TRUE), 
                function(x) 
                    read_csv(x, 
                        col_names = c("vehicle_id", "trip_id", "ts", "prior_mse", "posterior_mse", #"sumwt", "varwt",
                                      "post_speed", "prior_speed_var", "posterior_speed_var", "dist_to_path", 
                                      "Neff", "resample", "n_resample", "bad_sample"),
                        col_types = "ccidddddddiii", progress = FALSE) %>%
                    mutate(ts = as.POSIXct(ts, origin = "1970-01-01"))
            )
        ) %>% mutate(sim = sim, n_particles = siminfo[1], gps_error = siminfo[2], system_noise = siminfo[3])
    })
    if (inherits(x, "try-error")) return(NULL)
    saveRDS(x, file.path("simulations", sim, "modeleval.rds"))
    x
}

get_sims <- function(...) {
    do.call(bind_rows, pbapply::pblapply(list(...), get_sim_files))
}

sims <- list.files("simulations", pattern = "sim_")
sims <- sims[sapply(sims, function(s) file.exists(file.path("simulations", s, "timings.csv")))]
sims <- do.call(get_sims, as.list(sims))
sims <- sims %>% filter(dist_to_path >= 0 & dist_to_path < 50)

pdf("~/Dropbox/modeleval.pdf", onefile = TRUE, width = 14, height = 10)

## Observation from path
ggplot(sims %>% filter(dist_to_path >= 0 & dist_to_path < 30)) + 
    geom_density(aes(x = dist_to_path), fill = "gray")

## Effective sample size
ggplot(sims %>% filter(Neff <= n_particles & Neff > 0 & bad_sample == 0)) + 
    geom_violin(aes(x = factor(system_noise), Neff, fill = factor(system_noise))) +
    facet_grid(n_particles~gps_error, scales="free_y")

# ## sum of weights
# ggplot(sims) + 
#     geom_violin(aes(x = factor(system_noise), sumwt, fill = factor(system_noise))) +
#     facet_grid(n_particles~gps_error)

# ggplot(sims) + 
#     geom_violin(aes(x = factor(system_noise), varwt, fill = factor(system_noise))) +
#     facet_grid(n_particles~gps_error)


# ggplot(sims %>% filter(Neff < n_particles & Neff >= 0 & system_noise < 0.1 & gps_error == 3)) +
#     geom_point(aes(x = Neff, y = n_resample, color = factor(system_noise))) + 
#     facet_wrap(~ n_particles, scales = "free_x")

## Number of bad samples
ggplot(sims %>% group_by(n_particles, gps_error, system_noise) %>%
        summarize(p_bad = mean(bad_sample))) + 
    geom_col(aes(as.factor(n_particles), y = p_bad, fill = as.factor(system_noise)), position = "dodge") +
    facet_grid( ~ gps_error) +
    xlab("Numer of particles") + ylab("Rate of degeneration (proportion of iterations)")

## Number of Neff is good
ggplot(sims %>% group_by(n_particles, gps_error, system_noise) %>%
        summarize(p_good = mean(Neff >= 1000))) + 
    geom_col(aes(as.factor(n_particles), y = p_good, fill = as.factor(system_noise)), position = "dodge") +
    facet_grid( ~ gps_error) +
    xlab("Numer of particles") + ylab("1 / Rexampling rate (Neff large enough)")


ggplot(sims %>% group_by(n_particles, gps_error, system_noise) %>%
        summarize(p_bad = mean(bad_sample), p_good = mean(Neff < 1000))) + 
    geom_point(aes(p_good, p_bad, shape = as.factor(gps_error), color = as.factor(system_noise))) +
    facet_wrap( ~ n_particles) +
    xlab("Resampling Rate") + ylab("Degeneration Rate")

## Number of consecutive resamples
ggplot(sims) + 
    geom_step(aes(ts, n_resample, group = vehicle_id, color = as.factor(n_particles))) +
    facet_grid(gps_error ~ system_noise)

ggplot(sims) +
    geom_violin(aes(x = factor(system_noise), n_resample, fill = factor(system_noise))) +
    facet_grid(gps_error  ~ n_particles)

# ggplot(sims %>% group_by(n_particles, gps_error, system_noise) %>%
#         summarize(p_resample = mean(resample))) + 
#     geom_col(aes(as.factor(n_particles), y = p_resample, fill = as.factor(system_noise)), position = "dodge") +
#     facet_wrap( ~ gps_error)


## Distance from sample to obs
ggplot(sims %>% filter(prior_mse < 1000 & prior_mse > 0)) + 
    geom_violin(aes(factor(n_particles), prior_mse, fill = factor(n_particles))) +
    facet_grid(gps_error~system_noise)

ggplot(sims %>% filter(posterior_mse < 200 & posterior_mse > 0 & bad_sample == 0 & n_resample > 1)) +
    geom_violin(aes(factor(n_particles), pmin(posterior_mse / dist_to_path, 5), fill = factor(n_particles))) +
    facet_grid(gps_error~system_noise)


ggplot(sims %>% filter(dist_to_path < 50) %>% group_by(n_particles, gps_error, system_noise) %>%
        summarize(p_bad = mean(posterior_mse / dist_to_path > 2), p_good = mean(Neff < 1000))) + 
    geom_point(aes(p_good, p_bad, shape = as.factor(gps_error), color = as.factor(system_noise))) +
    facet_wrap( ~ n_particles, scales = "free") +
    xlab("Resampling Rate") + ylab("Proportion Posterior Mean Dist Too Large")

ggplot(sims %>% filter(dist_to_path < 20) %>% group_by(n_particles, gps_error, system_noise) %>%
        summarize(p_bad = mean(posterior_mse / dist_to_path > 2), p_good = mean(bad_sample))) + 
    geom_point(aes(p_good, p_bad, shape = as.factor(gps_error), color = as.factor(system_noise))) +
    facet_wrap( ~ n_particles) +
    xlab("Degeneration Rate") + ylab("Proportion Posterior Mean Dist Too Large")


## speed?
ggplot(sims %>% filter(n_resample > 1 & bad_sample == 0)) +
    geom_violin(aes(factor(n_particles), post_speed, fill = factor(n_particles))) +
    facet_grid(gps_error~system_noise)

## Speed variance
ggplot(sims %>% filter(post_speed > 0 & prior_speed_var > 0 & n_particles == 1000 & n_resample > 1 & bad_sample == 0)) + 
    geom_point(aes(post_speed, prior_speed_var)) +
    facet_grid(gps_error~system_noise)

ggplot(sims %>% filter(post_speed > 0 & prior_speed_var > 0 & n_particles == 1000 & n_resample > 1 & bad_sample == 0)) + 
    geom_violin(aes(factor(n_particles), sqrt(prior_speed_var), fill = factor(n_particles))) +
    facet_grid(gps_error~system_noise)

ggplot(sims %>% filter(post_speed > 0 & posterior_speed_var > 0 & n_particles == 1000 & n_resample > 1 & bad_sample == 0)) + 
    geom_point(aes(post_speed, prior_speed_var), col = 'gray') +
    geom_point(aes(post_speed, posterior_speed_var)) +
    facet_grid(gps_error~system_noise)

ggplot(sims %>% filter(posterior_speed_var > 0 & n_resample > 1 & bad_sample == 0 & posterior_speed_var < 10 & prior_speed_var < 10)) + 
    geom_violin(aes(factor(n_particles), (prior_speed_var)), color = 'gray') +
    geom_violin(aes(factor(n_particles), (posterior_speed_var), fill = factor(n_particles))) +
    facet_grid(gps_error~system_noise)


ggplot(sims %>% filter(gps_error == 3 & posterior_speed_var < 50 & n_resample > 1 & bad_sample == 0)) +
    geom_point(aes(Neff, posterior_speed_var), alpha = 0.1) +
    facet_grid(system_noise~n_particles, scales = "free_x")


dev.off()

#ggplot(sims %>% filter(posterior_mse < 200 & dist_to_path < 100)) +
#    geom_point(aes(dist_to_path, posterior_mse / prior_mse)) +
#    facet_grid(gps_error~system_noise)
    
#ggplot(sims %>% filter(prior_mse < 1000 & posterior_mse < 200)) +
#    geom_point(aes(prior_mse, posterior_mse)) +
#    facet_grid(sim~.)

# ggplot(sims %>% filter(posterior_mse < 100)) +
#     geom_hex(aes(ts, posterior_mse)) +
#     facet_grid(sim~.)

#ggplot(sims %>% filter(prior_mse < 1000 & posterior_mse < 200 & dist_to_path < 100)) +
#    geom_violin(aes(x = factor(sim), y = posterior_mse / prior_mse))



get_nw_times <- function(sim) {
    siminfo <- strsplit(sim, "_")[[1]][-1]
    if (grepl("e", siminfo[3])) siminfo[3] <- format(as.numeric(siminfo[3]), scientific = FALSE)
    siminfo <- as.numeric(gsub("-", ".", siminfo))
    read_csv(file.path("simulations", sim, "segment_states.csv"),
        col_names = c("segment_id", "timestamp", "mean", "var"), col_types = "cinn") %>%
        mutate(timestamp = as.POSIXct(timestamp, origin = "1970-01-01"),
               segment_id = factor(segment_id),
               sim = sim, n_particles = siminfo[1], gps_error = siminfo[2], system_noise = siminfo[3])
}

get_all_nw_times <- function() {
    sims <- list.files("simulations", pattern = "sim_")
    sims <- sims[sapply(sims, function(s) file.exists(file.path("simulations", s, "timings.csv")))]
    lapply(sims, get_nw_times) %>% bind_rows
}

nwtimes <- get_all_nw_times()

## one segment



# library(lme4)
# fit <- lmer(var ~ (n_particles + gps_error + system_noise | segment_id), data = nwtimes)

nwtimes.smry1 <- nwtimes %>% group_by(n_particles, gps_error, system_noise, segment_id) %>%
    summarize(sd.time = sd(var, na.rm = TRUE), mean.var = mean(var, na.rm = TRUE)) %>% ungroup

nwtimes.smry <- nwtimes.smry1 %>% group_by(n_particles, gps_error, system_noise) %>%
    summarize(mean = mean(sd.time, na.rm = TRUE), sd = mean(mean.var, na.rm = TRUE)) %>% ungroup()

# ggplot(nwtimes.smry) +
#     geom_point(aes(sqrt(mean.var/n_particles), var.var, color = as.factor(system_noise))) +
#     facet_grid(~n_particles)

ggplot(nwtimes.smry) +
    geom_point(aes(mean / sqrt(n_particles), sd / sqrt(n_particles), 
        shape = as.factor(system_noise))) +
    facet_grid(gps_error~n_particles) +
    xlab("Standard deviation of travel times") + ylab("Standard deviation of travel time uncertainty") +
    labs(color = "GPS Error (m)", shape = "System noise")

ggplot(nwtimes.smry) +
    geom_point(aes(as.numeric(as.factor(gps_error)) + as.numeric(as.factor(system_noise)) / 5, 
        mean.var/n_particles, color = as.factor(system_noise))) +
    facet_grid(~n_particles+gps_error)

ggplot(nwtimes.smry) +
    geom_point(aes(as.factor(system_noise), mean.var, color = as.factor(gps_error))) +
    facet_grid(~n_particles)

ggplot(nwtimes.smry) +
    geom_point(aes(gps_error, sqrt(mean.var/n_particles))) +
    facet_grid(~n_particles)





myacf <- function(x) {
    # x <- x[-(1:2)]
    cor(x[-1], x[-length(x)])
}
nwtimes.smry <- nwtimes %>% group_by(n_particles, gps_error, system_noise, segment_id) %>%
    summarize(autocov = myacf(mean)) %>% ungroup

    # summarize(speed.mean = mean(var), speed.var = var(var))

ggplot(nwtimes.smry %>% group_by(gps_error, n_particles, system_noise) %>%
        summarize(autocov = mean(autocov))) + 
    geom_point(aes(as.factor(n_particles), autocov, color = as.factor(gps_error))) +
    facet_grid( ~ system_noise)

# lm(autocov ~ n_particles + poly(system_noise, 2) + poly(gps_error, 2), data = nwtimes.smry) %>% summary



ggplot(nwtimes) +
    geom_point(aes(mean, var))

ggplot(nwtimes %>% filter(segment_id == "3073")) +
    geom_path(aes(timestamp, mean, color = as.factor(gps_error))) +
    facet_grid(as.factor(n_particles) ~ system_noise)


ggplot(nwtimes %>% filter(segment_id == "3073")) +
    geom_path(aes(timestamp, mean, color = as.factor(n_particles))) +
    facet_grid(system_noise ~ gps_error)





### OK OK OK so, actually, what is the RMSE for ETAs based on our super basic bitchin' algorithm?
source('scripts/common.R')

get_etas <- function(sim) {
    siminfo <- strsplit(sim, "_")[[1]][-1]
    if (grepl("e", siminfo[3])) siminfo[3] <- format(as.numeric(siminfo[3]), scientific = FALSE)
    siminfo <- as.numeric(gsub("-", ".", siminfo))
    all_sims(sim) %>%
        mutate(sim = sim, n_particles = siminfo[1], gps_error = siminfo[2], system_noise = siminfo[3])
}

get_all_etas <- function() {
    sims <- list.files("simulations", pattern = "sim_")
    sims <- sims[sapply(sims, function(s) file.exists(file.path("simulations", s, "timings.csv")))]
    do.call(get_sims, as.list(sims))
}

# etas <- get_all_etas()

## just one simulation

etas <- loadsim("sim000", gsub("etas_|.pb", "", tail(list.files("simulations/sim000/etas", pattern = ".pb"), 1))) %>%
    filter(time > timestamp) %>%
    mutate(
        timestamp = as.integer(timestamp),
        time = as.integer(time),
        eta = time - timestamp,
        lower = q0.025 - timestamp,
        upper = q0.975 - timestamp,
    ) %>%
    select(vehicle_id, trip_id, route_id, timestamp, stop_sequence, time, eta, lower, upper) %>%
    filter(!is.na(time) & lower > 0 & upper > 0) %>%
    filter(abs(time - timestamp) < 60*60)

ggplot(etas, aes(eta/60, stop_sequence)) + geom_point()

## one trip ...
trip <- etas %>% filter(trip_id == unique(etas$trip_id)[5])
ggplot(trip, aes(eta/60, stop_sequence)) + 
    geom_segment(aes(x = lower/60, xend = upper/60, yend = stop_sequence)) +
    geom_point()


load("simulations/arrivaldata.rda")
etas <- get_etas("sim000")


## Baseline - schedule RMSE ~ time-until-arrival
t1 <- arrivaldata %>% filter(trip_id == arrivaldata$trip_id[1]) %>% arrange(stop_sequence, type)
sched <- get_schedule(t1$trip_id[1]) %>%
    gather(key = "schedule_type", value = "schedule_time", arrival_time, departure_time) %>%
    mutate(schedule_type = gsub("_time", "", schedule_type)) %>%
    arrange(stop_sequence, schedule_type)

as.time <- function(x, date) as.POSIXct(paste(date, x))
t1delay <- t1 %>% select(trip_id, stop_sequence, type, time) %>% 
    left_join(sched, by = c("trip_id", "stop_sequence", "type" = "schedule_type")) %>% 
    mutate(time = as.POSIXct(time), 
           schedule_time = as.time(schedule_time, format(time[1], "%Y-%m-%d")),
           delay = as.integer(time - schedule_time))

schedETAs <-
    as_tibble(
        expand.grid(stop_sequence = 2:max(sched$stop_sequence),
                    time = t1delay$time[-1])) %>%
        left_join(t1delay %>% ungroup %>% select(time, delay), by = "time") %>%
        left_join(t1delay %>% ungroup %>% select(stop_sequence, schedule_time), by = "stop_sequence") %>%
        mutate(eta = schedule_time + delay)


ggplot(schedETAs) + 
    geom_point(aes(time, delay, colour = as.factor(stop_sequence)))


arrdat <- arrivaldata %>% filter(type == "arrival") %>%
    mutate(id = paste(trip_id, stop_sequence, sep = ":")) %>%
    group_by(id) %>% do((.) %>% tail(1)) %>% ungroup()

etasc <- etas %>% mutate(id = paste(trip_id, stop_sequence, sep = ":")) %>%
    inner_join(arrdat %>% select(id, time), by = "id", suffix = c("_estimated", "_actual"))

    #  %>%
    # inner_join(schedETAs %>% select())

ed <- etasc %>% filter(!is.na(time_actual) & !is.na(time_estimated)) %>%
    mutate(arrives_in = as.integer(time_actual - timestamp)) %>%
    filter(arrives_in > 0 & arrives_in < 60*60) %>%
    mutate(
        eta = as.integer(time_estimated - timestamp),
        eta_error = eta - arrives_in
    )
 # %>% filter(!is.na(q50) & !is.na(time_actual)) %>%
 #    filter(time_actual >= timestamp & abs(q50 - as.numeric(time_actual)) < 60*60) %>%
 #    mutate(in_pred_int = time_actual > q0 & time_actual < q100)

ggplot(ed %>% filter(abs(eta_error) < 2*60*60 & stop_sequence < 10)) +
    geom_histogram(aes(eta_error/60)) +
    facet_wrap(~stop_sequence)
    

ggplot(ed %>% filter(time_actual - timestamp < 60*60)) + geom_bar(aes(in_pred_int))

ggplot(ed) +
    geom_point(aes(as.integer(time_actual - timestamp)/60,))

ggplot(ed) +
    geom_point(aes(timestamp, (q50 - as.numeric(time_actual)) / 60), colour = 'gray') +
    geom_point(aes(timestamp, (q50 - as.numeric(time_actual)) / 60), data = ed %>% filter(in_pred_int)) +
    xlab("Time") +
    ylab("ETA error (minutes)")

## in 95% CI
ggplot(ed) +
    geom_point(aes(as.integer(time_actual - timestamp)/60, 
        ifelse(in_pred_int, 0, ifelse(time_actual < q0, as.numeric(time_actual) - q0,
                                        as.numeric(time_actual) - q100)) / 60)) +
        # (as.numeric(time_actual) - q5) / 60)) +
    ylab("ETA error (minutes)") +
    xlab("Time until arrival (min)")
    # geom_point(aes(timestamp, as.integer(time_actual - timestamp) / 60))

####### Timings
