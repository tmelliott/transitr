library(tidyverse)

get_sim_files <- function(sim) {
    x <- try({
        siminfo <- as.numeric(gsub("-", ".", strsplit(sim, "_")[[1]][-1]))
        do.call(bind_rows, 
            lapply(list.files(file.path("simulations", sim, "modeleval"), pattern="vehicle_.*\\.csv", full.names = TRUE), 
                function(x) 
                    read_csv(x, 
                        col_names = c("vehicle_id", "trip_id", "ts", "prior_mse", "posterior_mse", 
                                      "prior_speed_var", "posterior_speed_var", "dist_to_path", "Neff", "resample"),
                        col_types = "cciddddddi", progress = FALSE) %>%
                    mutate(ts = as.POSIXct(ts, origin = "1970-01-01"))
            )
        ) %>% mutate(sim = sim, n_particles = siminfo[1], gps_error = siminfo[2], system_noise = siminfo[3])
    })
    if (inherits(x, "try-error")) return(NULL)
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
ggplot(sims %>% filter(Neff <= n_particles & Neff >= 0)) + 
    geom_violin(aes(x = factor(system_noise), Neff, fill = factor(system_noise))) +
    facet_grid(n_particles~gps_error, scales="free_y")

## Distance from sample to obs
ggplot(sims %>% filter(dist_to_path < 200 & prior_mse < 5000)) + 
    geom_violin(aes(factor(n_particles), prior_mse, fill = factor(n_particles)), alpha = 0.5) +
    facet_grid(gps_error~system_noise)

ggplot(sims %>% filter(dist_to_path < 200) %>% mutate(posterior_mse = pmin(100, posterior_mse))) +
    geom_violin(aes(factor(n_particles), posterior_mse, fill = factor(n_particles)), alpha = 0.5) +
    facet_grid(system_noise~gps_error)

## Speed variance
ggplot(sims) + 
    geom_violin(aes(factor(n_particles), prior_speed_var, fill = factor(n_particles)), alpha = 0.5) +
    facet_grid(gps_error~system_noise)

ggplot(sims) + 
    geom_violin(aes(factor(n_particles), posterior_speed_var, fill = factor(n_particles)), alpha = 0.5) +
    facet_grid(gps_error~system_noise)


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




### OK OK OK so, actually, what is the RMSE for ETAs based on our super basic bitchin' algorithm?
