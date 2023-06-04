context("Create and Load GTFS database")
library(RSQLite)

tmpdb <- NULL
test_that("database created successfully", {
    nw <- create_gtfs()
    tmpdb <<- nw$database
    expect_is(nw, "trgtfs")
})

test_that("database loaded successfully", {
    nw <- load_gtfs(tmpdb)
    expect_is(nw, "trgtfs")
    expect_error(load_gtfs(tempfile()))
})


test_that("duplicate tables are detected", {
    expect_error(create_agency(tmpdb))
    expect_error(create_routes(tmpdb))
    expect_error(create_trips(tmpdb))
    expect_error(create_shapes(tmpdb))
    expect_error(create_stops(tmpdb))
    expect_error(create_stop_times(tmpdb))
    expect_error(create_calendar(tmpdb))
    expect_error(create_calendar_dates(tmpdb))
    expect_error(create_versions(tmpdb))
})


test_that("version API works", {
    if (Sys.getenv('APIKEY') == "") skip("No API key")

    nw <- load_gtfs(tmpdb) %>%
        version_api("https://api.at.govt.nz/v2/gtfs/versions")
    expect_warning(update_versions(nw))
    
    nw <- load_gtfs(tmpdb) %>%
        version_api("https://api.at.govt.nz/v2/gtfs/versions",
                    with_headers('Ocp-Apim-Subscription-Key' = Sys.getenv('APIKEY')))
    expect_true(has_version_api(nw))
    expect_output(update_versions(nw), "Versions updated")

    expect_error(send(1))
})


url <- "https://github.com/tmelliott/transitr/raw/develop/inst/extdata/auckland_gtfs.zip"
fzip <- tempfile(fileext = ".zip")
fdir <- file.path(tempdir(), "data")

skip <- inherits(try(download.file(url, fzip, quiet = TRUE)), "try-error")
if (!skip) unzip(fzip, exdir = fdir)

test_that("database updates from directory", {
    if (skip) skip("Couldn't download file")
    expect_is(create_gtfs(fdir, quiet = TRUE), "trgtfs")
    expect_output(create_gtfs(fdir))
})

test_that("database updates from a zip file", {
    if (skip) skip("Couldn't download file")
    expect_is(create_gtfs(fzip, quiet = TRUE), "trgtfs")
})

test_that("database updates from a URL (remote ZIP file)", {
    x <- try(create_gtfs(url, quiet = TRUE))
    if (inherits(x, "try-error")) skip("Couldn't download file")
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
