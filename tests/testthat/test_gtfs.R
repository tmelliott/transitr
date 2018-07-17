context("Create and Load GTFS database")
library(RSQLite)


tmpdb <- NULL
test_that("database created successfully", {
    nw <- create_gtfs()
    tmpdb <<- attr(nw$connection, "dbname")
    expect_is(nw, "trgtfs")
    dbDisconnect(nw$connection)
})

test_that("database loaded successfully", {
    con <- RSQLite::dbConnect(SQLite(), tmpdb)
    nw <- load_gtfs(con)
    expect_is(nw, "trgtfs")
    RSQLite::dbDisconnect(con)
    con <- RSQLite::dbConnect(SQLite(), tempfile())
    expect_error(load_gtfs(con))
    RSQLite::dbDisconnect(con)
})

test_that("duplicate tables are detected", {
    con <- RSQLite::dbConnect(SQLite(), tmpdb)
    expect_error(create_agency(con))
    expect_error(create_routes(con))
    expect_error(create_trips(con))
    expect_error(create_shapes(con))
    expect_error(create_stops(con))
    expect_error(create_stop_times(con))
    expect_error(create_calendar(con))
    expect_error(create_calendar_dates(con))
    expect_error(create_versions(con))
    RSQLite::dbDisconnect(con)
})


test_that("version API works", {
    if (Sys.getenv('APIKEY') == "") skip()
    con <- RSQLite::dbConnect(SQLite(), tmpdb)

    nw <- load_gtfs(con) %>%
        version_api("https://api.at.govt.nz/v2/gtfs/versions")
    expect_warning(update_versions(nw))
    
    nw <- load_gtfs(con) %>%
        version_api("https://api.at.govt.nz/v2/gtfs/versions",
                    headers =
                        list('Ocp-Apim-Subscription-Key' = Sys.getenv('APIKEY')))
    expect_true(has_version_api(nw))
    expect_output(update_versions(nw), "Versions updated")

    RSQLite::dbDisconnect(con)

    expect_error(send(1))
})


url <- "https://github.com/tmelliott/transitr/raw/develop/inst/extdata/auckland_gtfs.zip"
fzip <- tempfile(fileext = ".zip")
fdir <- file.path(tempdir(), "data")

download.file(url, fzip, quiet = TRUE)
unzip(fzip, exdir = fdir)

test_that("database updates from directory", {
    expect_is(create_gtfs(fdir, quiet = TRUE), "trgtfs")
    expect_output(create_gtfs(fdir))
})

test_that("database updates from a zip file", {
    expect_is(create_gtfs(fzip, quiet = TRUE), "trgtfs")
})

test_that("database updates from a URL (remote ZIP file)", {
    expect_is(create_gtfs(url, quiet = TRUE), "trgtfs")
})

test_that("table names are dealt with", {
    expect_equal(get_table_name("path/to/trips.txt"), "trips")
    expect_equal(get_table_name("path/to/routes.csv"), "routes")
    expect_warning(get_table_name("path/to/shapes.xlsx"))
    expect_true(check_valid_table("trips"))
    expect_true(check_valid_table("stop_times"))
    expect_false(check_valid_table("hello"))
})
