ca <- commandArgs(trailingOnly = TRUE)

if (length(ca) == 0)
    stop("You need to specify which simulation to run\n  `--args simN`\n")

library(transitr)
library(magrittr)

simdir <- file.path("simulations", ca[1])

setwd(simdir)

config <- jsonlite::read_json("config.json")
if (!is.null(config$simulation_history)) {
    if (dir.exists("history")) unlink("history", recursive = TRUE, force = TRUE)
    system(sprintf("rm -rf %s/*", config$simulation_history))
    system(sprintf("ln -sf %s history", config$simulation_history))
} else {
    if (dir.exists("history")) unlink("history", recursive = TRUE, force = TRUE)
    dir.create("history")
}

if (dir.exists("etas")) unlink("etas", recursive = TRUE, force = TRUE)
dir.create("etas")
if (file.exists("segment_states.csv")) unlink("segment_states.csv")

nw <- load_gtfs("../../fulldata.db", output = "etas.pb") %>%
    realtime_feed(sprintf("http://localhost:3000/%s/vehicle_positions", ca[1]), 
                  response = "protobuf")
nw <- do.call(set_parameters, c(list(nw), config))

## restart the server for this ID
RCurl::getURL(sprintf("localhost:3000/%s/reset", ca[1]))

model(nw)
