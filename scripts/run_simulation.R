ca <- commandArgs(trailingOnly = TRUE)

if (length(ca) == 0)
    stop("You need to specify which simulation to run\n  `--args simN`\n")

library(transitr)
library(magrittr)

simdir <- file.path("simulations", ca[1])

setwd(simdir)

nw <- load_gtfs("../../fulldata.db", output = "etas.pb") %>%
    realtime_feed(sprintf("http://localhost:3000/%s/vehicle_positions", ca[1]), 
                  response = "protobuf")
nw <- do.call(set_parameters, c(list(nw), jsonlite::read_json("config.json")))

## restart the server for this ID
RCurl::getURL(sprintf("localhost:3000/%s/reset", ca[1]))

model(nw)
