context("Set up real-time model")
library(RSQLite)

nw <- create_gtfs()

test_that("Parameters are set with valid defaults", {
    expect_is(nw$parameters, "list")

    expect_is(nw$parameters$n_core, "integer")

    expect_is(nw$parameters$n_particles, "integer")

    expect_is(nw$parameters$noise_model, "integer")

    expect_is(nw$parameters$system_noise, "numeric")
    expect_true(nw$parameters$system_noise > 0)

    expect_is(nw$parameters$pr_stop, "numeric")
    expect_true(nw$parameters$pr_stop >= 0 && nw$parameters$pr_stop <= 1)

    expect_is(nw$parameters$dwell_time, "numeric")
    expect_true(nw$parameters$dwell_time > 0)

    expect_is(nw$parameters$dwell_time_var, "numeric")
    expect_true(nw$parameters$dwell_time_var > 0)

    expect_is(nw$parameters$gamma, "numeric")
    expect_true(nw$parameters$gamma > 0)

    expect_is(nw$parameters$gps_error, "numeric")
    expect_true(nw$parameters$gps_error > 0)

    expect_is(nw$parameters$arrival_error, "numeric")
    expect_true(nw$parameters$arrival_error > 0)

    expect_is(nw$parameters$departure_error, "numeric")
    expect_true(nw$parameters$departure_error > 0)

    expect_is(nw$parameters$nw_system_noise, "numeric")
    expect_true(nw$parameters$nw_system_noise > 0)

    expect_is(nw$parameters$nw_measurement_error, "numeric")
    expect_true(nw$parameters$nw_measurement_error > 0)

    expect_is(nw$parameters$eta_model, "integer")

    expect_is(nw$parameters$save_timings, "logical")

    expect_is(nw$parameters$reset_method, "integer")
    expect_true(nw$parameters$reset_method > 0)
})

test_that("Parameters can be changed", {
    nw <- nw %>% set_parameters(
        n_core = 2L,
        n_particles = 20L,
        noise_model = 1L,
        system_noise = 1,
        pr_stop = 0.9,
        dwell_time = 15,
        dwell_time_var = 5,
        gamma = 2,
        gps_error = 20,
        arrival_error = 10,
        departure_error = 12,
        nw_system_noise = 0.002,
        nw_measurement_error = 50,
        eta_model = 1L,
        save_timings = TRUE,
        reset_method = 2L
    )

    expect_equal(
        nw$parameters,
        list(
            n_core = 2L,
            n_particles = 20L,
            noise_model = 1L,
            system_noise = 1,
            pr_stop = 0.9,
            dwell_time = 15,
            dwell_time_var = 5,
            gamma = 2,
            gps_error = 20,
            arrival_error = 10,
            departure_error = 12,
            nw_system_noise = 0.002,
            nw_measurement_error = 50,
            eta_model = 1L,
            save_timings = TRUE,
            reset_method = 2L
        )
    )
})
