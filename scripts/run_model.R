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
    realtime_feed("https://dl.dropboxusercontent.com/s/1fvto9ex649mkri/vehicle_locations.pb?dl=1", 
                  response = "protobuf")

model(nw, 500, ifelse(Sys.info()["nodename"] == "certellprd01", 6, 2))
