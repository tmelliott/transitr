library(transitr)

library(magrittr)

  if (file.exists("fulldata.db")) {
      nw <- load_gtfs("fulldata.db")
  } else {
      nw <- create_gtfs("https://cdn01.at.govt.nz/data/gtfs.zip",
                        db = "fulldata.db") %>%
      construct()     
  }

nw <- nw %>%
    realtime_feed("https://api.at.govt.nz/v2/public/realtime/vehiclelocations",
                  with_headers('Ocp-Apim-Subscription-Key' = Sys.getenv('APIKEY')),
                  response = "protobuf")

model(nw, 5000, 2)
