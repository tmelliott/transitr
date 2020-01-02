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

if (dir.exists("modeleval")) unlink("modeleval", recursive = TRUE, force = TRUE)
dir.create("modeleval")

if (file.exists("segment_states.csv")) unlink("segment_states.csv")
if (file.exists("segment_observations.csv")) unlink("segment_observations.csv")
if (file.exists("particle_travel_times.csv")) unlink("particle_travel_times.csv")
if (file.exists("eta_state.csv")) unlink("eta_state.csv")
if (file.exists("eta_quantiles.csv")) unlink("eta_quantiles.csv")
if (file.exists("vehicle_states.csv")) unlink("vehicle_states.csv")
if (file.exists("trip_vehicle_states.csv")) unlink("trip_vehicle_states.csv")

nw <- load_gtfs("../../at_gtfs.db", output = "etas.pb") %>%
    realtime_feed(c(sprintf("http://localhost:3000/%s/vehicle_positions", ca[1]),
    # realtime_feed(c(sprintf("http://localhost:3000/%s/50/100/vehicle_positions", ca[1]),
    # realtime_feed(c(sprintf("http://localhost:3000/%s/1566165600/minutes/0/vehicle_positions", ca[1]),
                    sprintf("http://localhost:3000/%s/trip_updates", ca[1])),
                  response = "protobuf")
nw <- do.call(set_parameters, c(list(nw), config))

## restart the server for this ID
RCurl::getURL(sprintf("localhost:3000/%s/reset", ca[1]))

model(nw)
