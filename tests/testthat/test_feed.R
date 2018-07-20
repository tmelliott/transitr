context("Connect to realtime GTFS feed")

nw <- create_gtfs(system.file("extdata", "auckland_gtfs.zip", package = "transitr"),
                  quiet = TRUE) %>%
    construct()

test_that("URL gets added to object", {
    if (Sys.getenv('APIKEY') == "") skip("No API key")
    n1 <- nw %>% realtime_feed("https://api.at.govt.nz/v2/public/realtime/vehiclelocations",
                         with_headers('Ocp-Apim-Subscription-Key' = Sys.getenv('APIKEY')))   
    n2 <- nw %>% realtime_feed("https://api.at.govt.nz/v2/public/realtime/vehiclelocations",
                         with_headers('Ocp-Apim-Subscription-Key' = Sys.getenv('APIKEY')),
                         response = "protobuf")
    expect_is(n1$apis$realtime, "trapi")
})