context("Construct the GTFS road network")
library(RSQLite)

nw <- create_gtfs(system.file("extdata", "auckland_gtfs.zip", package = "transitr"),
                  quiet = TRUE)

test_that("network shapes makes sense", {
    sh <- load_shapes(nw)
    expect_is(sh, "network.shape.list")
    expect_equal(length(sh), 17)
})

test_that("network gets constructed correctly", {
    n <- construct(nw)
    expect_equal(dim(load_road_segments(nw)), c(379, 4))
    expect_equal(dim(load_intersections(nw)), c(369, 3))
    expect_equal(dim(load_shape_segments(nw)), c(581, 4))
})
