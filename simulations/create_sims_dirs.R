config <- list(
    n_core = 4,
    n_particles = 2000,
    noise_model = 2,
    system_noise = c(0.6, 0.7, 0.8, 0.9),
    pr_stop = 0.1,
    dwell_time = 6.0,
    gamma = 5.0,
    gps_error = 3,
    arrival_error = 5,
    departure_error = 5,
    nw_system_noise = 0.001,
    nw_measurement_error = 500,
    save_timings = TRUE
    )

grid <- do.call(expand.grid, config)

for (i in 1:nrow(grid)) {
    conf <- grid[i, , drop = TRUE]
    dirname <- glue::glue("sim_3341-{conf$system_noise}")
    if (!dir.exists(dirname))
        dir.create(dirname)
    if (!file.exists(file.path(dirname, "config.json")))
        jsonlite::write_json(conf, 
            file.path(dirname, "config.json"), 
            auto_unbox = TRUE, 
            pretty = TRUE
        )
}

# existing <- list.files(pattern = 'sim_', include.dirs = TRUE)
# todelete <- existing[!existing %in% dirs]
# sapply(todelete, unlink, recursive = TRUE, force = TRUE)


