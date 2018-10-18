library(tidyverse)

get_sim_files <- function(sim) {
    do.call(bind_rows, 
        pbapply::pblapply(list.files(file.path("simulations", sim, "modeleval"), pattern="vehicle_.*\\.csv", full.names = TRUE), 
            function(x) 
                read_csv(x, 
                    col_names = c("vehicle_id", "trip_id", "ts", "prior_mse", "posterior_mse", "dist_to_path", "Neff", "resample"),
                    col_types = "cciddddi", progress = FALSE) %>%
                mutate(ts = as.POSIXct(ts, origin = "1970-01-01"))
        )
    ) %>% mutate(sim = sim)
}

get_sims <- function(...) {
    do.call(bind_rows, lapply(list(...), get_sim_files))
}

sims <- get_sims("sim000", "sim001", "sim002", "sim003", "sim004")

## Effective sample size
ggplot(sims %>% filter(Neff <= 3000 & Neff >= 0)) + 
    geom_histogram(aes(x = Neff), bins = 100) +
    facet_grid(sim~.)


## Observation from path
ggplot(sims %>% filter(dist_to_path >= 0 & dist_to_path < 30 & sim == "sim000")) + 
    geom_density(aes(x = dist_to_path), fill = "gray")

## Distance from sample to obs
ggplot(sims %>% mutate(posterior_mse = pmin(200, posterior_mse))) +
    geom_histogram(aes(x = posterior_mse), bins = 50) +
    facet_grid(sim~.)

ggplot(sims %>% filter(posterior_mse < 200)) +
    geom_point(aes(dist_to_path, posterior_mse / prior_mse)) +
    facet_grid(sim~.)
    # geom_abline(colour = 'orangered') +
    
ggplot(sims %>% filter(prior_mse < 1000 & posterior_mse < 200)) +
    geom_point(aes(prior_mse, posterior_mse)) +
    facet_grid(sim~.)

# ggplot(sims %>% filter(posterior_mse < 100)) +
#     geom_hex(aes(ts, posterior_mse)) +
#     facet_grid(sim~.)

ggplot(sims %>% filter(prior_mse < 1000 & posterior_mse < 200 & dist_to_path < 100)) +
    geom_violin(aes(x = factor(sim), y = posterior_mse / prior_mse))
