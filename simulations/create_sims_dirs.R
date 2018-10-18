n <- c(1000, 2000, 5000, 10000)
gps <- c(3, 5, 8, 10)
noise <- c(0.01, 0.02, 0.05, 0.1, 0.5)

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
apply(grid, 1, function(x) {
    config$n_particles <- x[["n"]]
    config$gps_error <- x[["gps"]]
    config$system_noise <- x[["noise"]]
    dirname <- glue::glue("sim_{config$n_particles}_{config$gps_error}_{config$system_noise}")
    dirname <- gsub("\\.", "-", dirname)
    if (dir.exists(dirname)) unlink(dirname, TRUE, TRUE)
    dir.create(dirname)
    jsonlite::write_json(config, file.path(dirname, "config.json"), auto_unbox = TRUE, pretty = TRUE)
})
