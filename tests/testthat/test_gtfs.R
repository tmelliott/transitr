context("Create and Load GTFS database")
library(RSQLite)

test_that("database created successfully", {
    nw <- create_gtfs()
    expect_is(nw, "trgtfs")

    dbDisconnect(nw$connection)
})

test_that("database loaded successfully", {
    con <- RSQLite::dbConnect(SQLite(), tempfile())
    nw <- load_gtfs(con)
    expect_is(nw, "trgtfs")
    RSQLite::dbDisconnect(con)
})


url <- "https://cdn01.at.govt.nz/data/gtfs.zip"
fzip <- tempfile(fileext = ".zip")
fdir <- file.path(tempdir(), "data")

download.file(url, fzip)
unzip(fzip, exdir = fdir)

nw <- create_gtfs(fdir)
test_that("database updates from directory", {
    expect_is(nw, "trgtfs")
})

test_that("database updates from a zip file", {
    expect_equal(nw, create_gtfs(fzip))
})

test_that("database updates from a URL (remote ZIP file)", {
    expect_equal(nw, create_gtfs(url))
})
