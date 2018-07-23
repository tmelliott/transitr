library(transitr)

library(magrittr)

nw <- create_gtfs(system.file("extdata", "auckland_gtfs.zip", 
                              package = "transitr"),
                  quiet = TRUE) %>%
    construct() %>%
    realtime_feed("https://api.at.govt.nz/v2/public/realtime/vehiclelocations",
                  with_headers('Ocp-Apim-Subscription-Key' = Sys.getenv('APIKEY')),
                  response = "protobuf")

model(nw)
