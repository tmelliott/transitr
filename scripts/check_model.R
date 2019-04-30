source("scripts/common.R")

sim <- "sim000"
vids <- unique(gsub(
    "^v|_mutate\\.csv", 
    "", 
    list.files(
        file.path("simulations", sim, "history"),
        pattern = "^v[A-Z0-9]*_mutate\\.csv"
    )
))

v_hist <- function(id, type = c("mutate", "update"), sim = "sim000") {
    type <- match.arg(type)
    f <- sprintf(
        "simulations/%s/history/v%s_%s.csv",
        sim, id, type
    )
    if (type == "mutate") {
        cn <- c("timestamp", "trip_id", "type", "obs_latitude", "obs_longitude",
            "stop_index", "dist_to_route", "delta", "distance", "speed",
            "latitude", "longitude", "dist_to_obs", "sum_llh")
        ct <- "iccnnininnnnnn"
    } else {
        cn <- c("timestamp", "trip_id", "type", "distance", "speed", "latitude", 
            "longitude", "dist_to_obs", "Neff", "keep")
        ct <- "iccnnnnnni"
    }
    read_csv(f, col_names = cn, col_types = ct)
}

vi <- 2
vv <- v_hist(vids[vi]) %>% 
    left_join(
        v_hist(vids[vi], "update") %>% 
            select(timestamp, trip_id, type, distance, speed, latitude, 
                longitude, dist_to_obs, Neff, keep),
        by = c("timestamp", "trip_id", "type"),
        suffix = c("", "_update")
    ) %>%
    mutate(time = as.POSIXct(timestamp, origin = "1970-01-01"))


ggplot(vv, aes(time)) +
    geom_point(aes(y = distance_update, colour = type, size = sum_llh)) +
    facet_wrap(~trip_id, scales = "free")

ggplot(vv, aes(time)) +
    geom_point(aes(y = dist_to_obs_update, colour = type, size = sum_llh)) +
    facet_wrap(~trip_id, scales = "free")

ggplot(vv %>% filter(type == "gps"), aes(time)) +
    geom_point(aes(y = dist_to_route)) +
    geom_point(aes(y = dist_to_obs_update, size = sum_llh), 
        colour = "orangered", alpha = 0.5) +
    facet_wrap(~trip_id, scales = "free")

sh <- do.call(
    bind_rows,
    lapply(unique(vv$trip_id), function(t) 
        getshape(trip = t) %>% mutate(trip_id = t))
)

ggplot(vv, aes(obs_longitude, obs_latitude)) +
    geom_path(aes(shape_pt_lon, shape_pt_lat, group = shape_id), data = sh) +
    geom_point() +
    geom_segment(aes(xend = longitude, yend = latitude), colour = "orangered") +
    geom_point(aes(longitude, latitude, colour = "orangered"), size = 1) +
    geom_segment(aes(xend = longitude_update, yend = latitude_update), colour = "blue") +
    geom_point(aes(longitude_update, latitude_update), colour = "blue", size = 1) +
    facet_wrap(~trip_id)




vi <- 2
vv2 <- v_hist(vids[vi], sim = "sim001") %>% 
    full_join(
        v_hist(vids[vi], "update", sim = "sim001") %>% 
            select(timestamp, trip_id, type, distance, speed, latitude, 
                longitude, dist_to_obs, Neff, keep),
        by = c("timestamp", "trip_id", "type"),
        suffix = c("", "_update")
    ) %>%
    mutate(time = as.POSIXct(timestamp, origin = "1970-01-01"))

vv2 %>% filter(type == "gps") %>% pluck("keep")

ggplot(vv2, aes(time)) +
    geom_point(aes(y = distance_update, colour = type, size = sum_llh)) +
    facet_wrap(~trip_id, scales = "free")

ggplot(vv2, aes(time)) +
    geom_point(aes(y = dist_to_obs_update, colour = type, size = sum_llh)) +
    facet_wrap(~trip_id, scales = "free")

ggplot(vv2 %>% filter(type == "gps"), aes(time)) +
    geom_point(aes(y = dist_to_route)) +
    geom_point(aes(y = dist_to_obs_update, size = sum_llh), 
        colour = "orangered", alpha = 0.5) +
    facet_wrap(~trip_id, scales = "free")

sh <- do.call(
    bind_rows,
    lapply(unique(vv2$trip_id), function(t) 
        getshape(trip = t) %>% mutate(trip_id = t))
)

ggplot(vv2, aes(obs_longitude, obs_latitude)) +
    geom_path(aes(shape_pt_lon, shape_pt_lat, group = shape_id), data = sh) +
    geom_point() +
    geom_segment(aes(xend = longitude, yend = latitude, alpha = delta), colour = "orangered") +
    geom_point(aes(longitude, latitude, colour = "orangered"), size = 1) +
    geom_segment(aes(xend = longitude_update, yend = latitude_update), colour = "blue") +
    geom_point(aes(longitude_update, latitude_update), colour = "blue", size = 1) +
    facet_wrap(~trip_id)


vv2 %>% group_by(timestamp) %>%
    summarize(
        n = n(),
        dmin = min(delta),
        dmax = max(delta)
    ) %>% filter(n>1)

ggplot(vv, aes(sum_llh)) + geom_histogram()
ggplot(vv2, aes(sum_llh)) + geom_histogram()


############## Just like the all of it
source("scripts/common.R")

results_skip <- 
    list.files(
        "simulations/sim000/history",
        pattern = "vehicle_[A-Z0-9]*.csv", 
        full.names = TRUE
    ) %>%
    lapply(
        read_csv,
        col_types = cols(
            vehicle_id = col_character(),
            event_timestamp = col_integer(),
            vehicle_timestamp = col_integer()
        )
    ) %>%
    bind_rows() %>%
    mutate(
        event_timestamp = as.POSIXct(event_timestamp, origin = "1970-01-01"),
        vehicle_timestamp = as.POSIXct(vehicle_timestamp, origin = "1970-01-01")
    )

results_revert <- 
    list.files(
        "simulations/sim001/history",
        pattern = "vehicle_[A-Z0-9]*.csv", 
        full.names = TRUE
    ) %>%
    lapply(
        read_csv,
        col_types = cols(
            vehicle_id = col_character(),
            event_timestamp = col_integer(),
            vehicle_timestamp = col_integer()
        )
    ) %>%
    bind_rows() %>%
    mutate(
        event_timestamp = as.POSIXct(event_timestamp, origin = "1970-01-01"),
        vehicle_timestamp = as.POSIXct(vehicle_timestamp, origin = "1970-01-01")
    )

mean(results_skip$state_type == "initialize")
mean(results_revert$state_type == "initialize")

sum(results_skip$action == "skip", na.rm = TRUE)
sum(results_revert$action == "revert_state", na.rm = TRUE)

sum(results_skip$action == "resample", na.rm = TRUE) / nrow(results_skip)
sum(results_revert$action == "resample", na.rm = TRUE) / nrow(results_revert)

results_skip %>% filter(state_type == "update" & event_type == "gps") %>%
    ggplot(aes(sum_llh)) + geom_histogram()
results_revert %>% filter(state_type == "update" & event_type == "gps") %>%
    ggplot(aes(sum_llh)) + geom_histogram()


## noise models

sims <- list.files("simulations", pattern = "sim_")

results <- pbapply::pblapply(sims,
    function(sim) {
        x <- try({
            list.files(
                sprintf("simulations/%s/history", sim),
                pattern = "vehicle_[A-Z0-9]*.csv", 
                full.names = TRUE
            ) %>%
            lapply(
                read_csv,
                col_types = cols(
                    vehicle_id = col_character(),
                    event_timestamp = col_integer(),
                    vehicle_timestamp = col_integer()
                )
            ) %>%
            bind_rows() %>%
            mutate(
                event_timestamp = as.POSIXct(event_timestamp, 
                    origin = "1970-01-01"),
                vehicle_timestamp = as.POSIXct(vehicle_timestamp, 
                    origin = "1970-01-01"),
                sim = sim
            )
        }, silent = TRUE)
        if (inherits(x, "try-error")) return(NULL)
        x
    }
) %>% bind_rows %>%
    mutate(
        event_type = as.factor(event_type),
        action = fct_explicit_na(action, 'none'),
        sim = as.factor(sim),
        noise_model = as.factor(ifelse(grepl("3340", sim), "1", "2")),
        noise = as.numeric(gsub(".+\\-", "", sim))
    )

# ggplot(results, aes(event_timestamp, sum_llh)) +
#     geom_point() +
#     facet_wrap(~sim, ncol = 1)

# ggplot(results %>% 
#         filter(vehicle_position_error < 100 &
#             state_type %in% c('initialize', 'mutate', 'update')), 
#     aes(event_timestamp, vehicle_position_error)) +
#     geom_point(aes(colour = state_type)) +
#     facet_wrap(~sim, ncol = 1)

# iNZightPlots::iNZightPlot(state_type, sim, data = results,
#     inference.type = "conf", g1 = noise_model)

# iNZightPlots::iNZightPlot(action, sim, data = results)

# iNZightPlots::iNZightPlot(Neff, sim, data = results)

smry <- results %>% group_by(noise_model, noise) %>%
    summarize(
        pr_init = mean(state_type == "initialize"),
        pr_mutate = mean(state_type == "mutate"),
        pr_revert = mean(state_type == "revert"),
        pr_update = mean(state_type == "update"),
        n = n()
    ) %>% 
    filter(noise > 0.4)

gridExtra::grid.arrange(
    ggplot(smry, aes(noise, colour = noise_model)) +
        facet_wrap(~noise_model, scales = "free") +
        geom_linerange(aes(
            ymin = pr_init - 1.96 * pr_init * (1 - pr_init) / sqrt(n),
            ymax = pr_init + 1.96 * pr_init * (1 - pr_init) / sqrt(n))
        ) +
        geom_point(aes(y = pr_init)),
    ggplot(smry, aes(noise, colour = noise_model)) +
        facet_wrap(~noise_model, scales = "free") +
        geom_linerange(aes(
            ymin = pr_mutate - 1.96 * pr_mutate * (1 - pr_mutate) / sqrt(n),
            ymax = pr_mutate + 1.96 * pr_mutate * (1 - pr_mutate) / sqrt(n))
        ) +
        geom_point(aes(y = pr_mutate)),
    ggplot(smry, aes(noise, colour = noise_model)) +
        facet_wrap(~noise_model, scales = "free") +
        geom_linerange(aes(
            ymin = pr_revert - 1.96 * pr_revert * (1 - pr_revert) / sqrt(n),
            ymax = pr_revert + 1.96 * pr_revert * (1 - pr_revert) / sqrt(n))
        ) +
        geom_point(aes(y = pr_revert)),
    ggplot(smry, aes(noise, colour = noise_model)) +
        facet_wrap(~noise_model, scales = "free") +
        geom_linerange(aes(
            ymin = pr_update - 1.96 * pr_update * (1 - pr_update) / sqrt(n),
            ymax = pr_update + 1.96 * pr_update * (1 - pr_update) / sqrt(n))
        ) +
        geom_point(aes(y = pr_update)),
    ncol = 1
)

## at the end of the day, its which method is better at
#  estimating travel time that wins
