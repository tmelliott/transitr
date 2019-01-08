context("Set up real-time model")
library(RSQLite)

nw <- create_gtfs()

test_that("Parameters are set with valid defaults", {
    expect_is(nw$parameters, "list")

    expect_is(nw$parameters$n_core, "integer")

    expect_is(nw$parameters$n_particles, "integer")

    expect_is(nw$parameters$system_noise, "numeric")
    expect_true(nw$parameters$system_noise > 0)
    
    expect_is(nw$parameters$pr_stop, "numeric")
    expect_true(nw$parameters$pr_stop >= 0 && nw$parameters$pr_stop <= 1)
    
    expect_is(nw$parameters$dwell_time, "numeric")
    expect_true(nw$parameters$dwell_time > 0)

    expect_is(nw$parameters$gamma, "numeric")
    expect_true(nw$parameters$gamma > 0)

    expect_is(nw$parameters$gps_error, "numeric")
    expect_true(nw$parameters$gps_error > 0)

    expect_is(nw$parameters$arrival_error, "numeric")
    expect_true(nw$parameters$arrival_error > 0)

    expect_is(nw$parameters$departure_error, "numeric")
    expect_true(nw$parameters$departure_error > 0)

    expect_is(nw$parameters$save_timings, "logical")
})

test_that("Parameters can be changed", {
    nw <- nw %>% set_parameters(
        n_core = 2L,
        n_particles = 20L,
        system_noise = 1,
        pr_stop = 0.9,
        dwell_time = 5,
        gamma = 2,
        gps_error = 20,
        arrival_error = 10,
        departure_error = 12,
        save_timings = TRUE
    )

    expect_equal(
        nw$parameters,
        list(
            n_core = 2L,
            n_particles = 20L,
            system_noise = 1,
            pr_stop = 0.9,
            dwell_time = 5,
            gamma = 2,
            gps_error = 20,
            arrival_error = 10,
            departure_error = 12,
            save_timings = TRUE
        )
    )
})