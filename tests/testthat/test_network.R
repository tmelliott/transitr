context("Construct the GTFS road network")
library(RSQLite)

nw <- create_gtfs(
    system.file("extdata", "auckland_gtfs.zip", package = "transitr"),
    quiet = TRUE
)

test_that("network shapes makes sense", {
    sh <- load_shapes(nw)
    expect_is(sh, "network.shape.list")
    expect_equal(length(sh), 15)
})

test_that("network tables get constructed", {
    create_network_tables(nw)
    con <- db_connect(nw$database)
    expect_true(dbExistsTable(con, "road_segments"))
    expect_true(dbExistsTable(con, "intersections"))
    expect_true(dbExistsTable(con, "shape_nodes"))
    expect_true(dbExistsTable(con, "nodes"))
})

test_that("network gets constructed correctly", {
    construct_network(nw, node_threshold = 0.0)

    expect_equal(dim(load_nodes(nw)), c(356, 4))
    expect_equal(
        sapply(load_nodes(nw), class),
        c(
            node_id = "integer", node_type = "integer",
            node_lon = "numeric", node_lat = "numeric"
        )
    )
    expect_equal(dim(load_road_segments(nw)), c(362, 4))
    expect_equal(dim(load_intersections(nw)), c(0, 3))
    expect_equal(dim(load_shape_nodes(nw)), c(569, 4))

    # node shapes should have same max as stops and shape
    expect_true(
        all(abs(with(load_shape_nodes(nw), tapply(distance_traveled, shape_id, max)) -
            sapply(load_shapes(nw), function(x) max(x[,3]))) < 0.05)
    )
})

## system(sprintf("cp %s tests/testthat/auckland_gtfs.db", nw$database))
