library(transitr)

library(magrittr)

if (Sys.getenv("USE_FULLDATA") != "") {
    if (file.exists("fulldata.db")) {
        nw <- load_gtfs("fulldata.db")
    } else {
        nw <- create_gtfs("https://cdn01.at.govt.nz/data/gtfs.zip",
                          db = "fulldata.db") %>%
        construct()     
    }
} else {
    nw <- create_gtfs(system.file("extdata", "auckland_gtfs.zip", 
                                  package = "transitr"), 
                      quiet = TRUE) %>%
        construct()    
}

nw <- nw %>%
    realtime_feed("https://api.at.govt.nz/v2/public/realtime/vehiclelocations",
                  with_headers('Ocp-Apim-Subscription-Key' = Sys.getenv('APIKEY')),
                  response = "protobuf")

model(nw, 5)
