library(transitr)

library(magrittr)

fulldata <- Sys.getenv("USE_FULLDATA") != ""
f <- ifelse(fulldata, "https://cdn01.at.govt.nz/data/gtfs.zip",
            system.file("extdata", "auckland_gtfs.zip", package = "transitr"))
nw <- create_gtfs(f, quiet = TRUE) %>%
    construct() %>%
    realtime_feed("https://api.at.govt.nz/v2/public/realtime/vehiclelocations",
                  with_headers('Ocp-Apim-Subscription-Key' = Sys.getenv('APIKEY')),
                  response = "protobuf")

model(nw)
