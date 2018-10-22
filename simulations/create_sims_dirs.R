n <- c(500, 1000, 2000, 4000, 8000)
gps <- c(2, 3, 5, 8)
noise <- c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1)

config <- list(
    n_core = 8,
    n_particles = NA,
    gps_error = NA,
    system_noise = NA,
    pr_stop = 0.1,
    dwell_time = 6.0,
    gamma = 5.0,
    save_timings = TRUE
    )

grid <- expand.grid(n = n, gps = gps, noise = noise)
dirs <- apply(grid, 1, function(x) {
    config$n_particles <- x[["n"]]
    config$gps_error <- x[["gps"]]
    config$system_noise <- x[["noise"]]
    dirname <- glue::glue("sim_{config$n_particles}_{config$gps_error}_{config$system_noise}")
    dirname <- gsub("\\.", "-", dirname)
    if (!dir.exists(dirname))
        dir.create(dirname)
    if (!file.exists(file.path(dirname, "config.json")))
        jsonlite::write_json(config, file.path(dirname, "config.json"), auto_unbox = TRUE, pretty = TRUE)
    dirname
})

existing <- list.files(pattern = 'sim_', include.dirs = TRUE)
todelete <- existing[!existing %in% dirs]
sapply(todelete, unlink, recursive = TRUE, force = TRUE)


