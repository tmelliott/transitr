ca <- commandArgs(trailingOnly = TRUE)

library(transitr)
library(magrittr)

simdir <- file.path("simulations", ca[1])

setwd(simdir)

nw <- load_gtfs("../../fulldata.db", output = "etas.pb") %>%
    realtime_feed("file:///home/emperor/Dropbox/gtfs/vehicle_locations.pb", 
                  response = "protobuf")
nw <- do.call(set_parameters, c(list(nw), jsonlite::read_json("config.json")))

model(nw)
