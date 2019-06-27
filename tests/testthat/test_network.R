context("Construct the GTFS road network")
library(RSQLite)

nw <- create_gtfs(
    system.file("extdata", "auckland_gtfs.zip", package = "transitr"),
    quiet = TRUE
)

test_that("network shapes makes sense", {
    sh <- load_shapes(nw)
    expect_is(sh, "network.shape.list")
    expect_equal(length(sh), 17)
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
    construct_network(nw)

    expect_equal(dim(load_nodes(nw)), c(369, 4))
    expect_equal(
        sapply(load_nodes(nw), class),
        c(
            node_id = "integer", node_type = "integer",
            node_lon = "numeric", node_lat = "numeric"
        )
    )
    expect_equal(dim(load_road_segments(nw)), c(379, 4))
    expect_equal(dim(load_intersections(nw)), c(0, 3))
    expect_equal(dim(load_shape_segments(nw)), c(598, 4))
})

## file.copy(nw$database, "auckland_gtfs.db")
