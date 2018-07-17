context("Construct the GTFS road network")
library(RSQLite)

nw <- create_gtfs(system.file("extdata", "auckland_gtfs.zip", package = "transitr"),
                  quiet = TRUE)

test_that("network gets constructed correctly", {
    expect_is(load_shapes(nw), "list")
})
