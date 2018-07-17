context("Construct the GTFS road network")
library(RSQLite)

nw <- create_gtfs(system.file("extdata", "auckland_gtfs.zip", package = "transitr"),
                  quiet = TRUE)

# test_that("network shapes makes sense", {
#     sh <- load_shapes(nw)
#     expect_is(sh, "network.shape.list")
#     expect_equal(length(sh), 17)
# })

# test_that("network gets constructed correctly", {
#     n <- construct(nw)
#     expect_is(n, "network")
# })

dbDisconnect(nw$connection)