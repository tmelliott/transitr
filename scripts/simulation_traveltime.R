source("scripts/common.R")

## --- load a route and its stops (for a trip)
db <- "at_gtfs.db"
con <- dbConnect(SQLite(), db)
route <- con %>% tbl("routes") %>% 
    filter(
        route_short_name == "NX1" & 
        route_long_name %like% "Albany Station%"
    ) %>% 
    arrange(desc(version)) %>%
    head(1) %>% collect()
trip <- con %>% tbl("trips") %>%
    filter(route_id == !!route$route_id) %>%
    head(1) %>% collect()
stops <- con %>% tbl("stop_times") %>%
    filter(trip_id == !!trip$trip_id) %>%
    select(
        trip_id, stop_id, stop_sequence, 
        arrival_time, departure_time, shape_dist_traveled
    ) %>%
    left_join(
        con %>% tbl("stops") %>% select(stop_id, stop_lat, stop_lon)
    ) %>%
    arrange(stop_sequence) %>% collect()
shape <- con %>% tbl("shapes") %>%
    filter(shape_id == !!trip$shape_id) %>%
    arrange(shape_pt_sequence) %>% collect() %>%
    mutate(
        shape_dist_traveled = c(0, 
            (cbind(shape_pt_lon, shape_pt_lat) %>%
                geosphere::distGeo())[-n()]
        ) %>% cumsum
    )
dbDisconnect(con)


## --- generate some simulated data with known travel times
make_protobuf_feed <- function(data) {
    new(transit_realtime.FeedMessage,
        header = new(transit_realtime.FeedHeader,
            gtfs_realtime_version = "1.0",
            incrementality = 0,
            timestamp = data$timestamp
        ),
        entity = list(
            new(transit_realtime.FeedEntity,
                id = "testentity",
                is_deleted = FALSE,
                # if data$type %in% c("arrival", "departure")
                trip_update = new(transit_realtime.TripUpdate,
                    trip = new(transit_realtime.TripDescriptor,
                        trip_id = data$trip_id,
                        route_id = data$route_id,
                        direction_id = 0,
                        start_time = stops$arrival_time[1],
                        start_date = Date,
                        schedule_relationship = 0
                    ),
                    vehicle = new(transit_realtime.VehicleDescriptor,
                        id = data$vehicle_id
                    ),
                    stop_time_update = list(
                        new(transit_realtime.TripUpdate.StopTimeUpdate,
                            stop_sequence = data$stop_sequence,
                            stop_id = data$stop_id,
                            ## if data$type == "departure"
                            arrival = if (data$type == "arrival") 
                                new(transit_realtime.TripUpdate.StopTimeEvent,
                                    time = data$time
                                ) else NULL,
                            departure = if (data$type == "departure")
                                new(transit_realtime.TripUpdate.StopTimeEvent,
                                    time = data$time
                                ) else NULL
                        )
                    ),
                    timestamp = data$time
                )
            )
        )
    )
}

Date <- "2019-04-01"
start <- as.integer(as.POSIXct(
    sprintf("%s %s", Date, stops$arrival_time[1])
))

speed <- 15 # constant speed of 15m/s for simulation 1
stops$shape_dist_traveled <- sapply(1:nrow(stops), function(i) {
    shape$shape_dist_traveled[
        which.min(geosphere::distGeo(
            stops[i, c("stop_lon", "stop_lat")] %>% as.matrix,
            shape[, c("shape_pt_lon", "shape_pt_lat")] %>% as.matrix
        ))
    ]
})

travel_times <- round(diff(stops$shape_dist_traveled) / speed)
dwell_times <- c(0, rep(60, length(travel_times) - 1), 0)
data <- rbind(
    data.frame(
        vehicle_id = "TEST",
        trip_id = trip$trip_id,
        route_id = route$route_id,
        timestamp = start + cumsum(c(0, travel_times[-length(travel_times)])) + 
            cumsum(dwell_times[-length(dwell_times)]),
        type = "departure",
        stop_sequence = stops$stop_sequence[-nrow(stops)],
        stop_id = stops$stop_id[-nrow(stops)],
        stringsAsFactors = FALSE
    ),
    data.frame(
        vehicle_id = "TEST",
        trip_id = trip$trip_id,
        route_id = route$route_id,
        timestamp = start + cumsum(travel_times) + cumsum(dwell_times[-length(dwell_times)]),
        type = "arrival",
        stop_sequence = stops$stop_sequence[-1],
        stop_id = stops$stop_id[-1],
        stringsAsFactors = FALSE
    )
) %>% arrange(timestamp)

simdir <- file.path("scripts", "simulation_traveltime_01")
dir <- "scripts/simulation_traveltime_01/archive"
unlink(simdir, TRUE, TRUE)
dir.create(dir, recursive = TRUE)
invisible(lapply(
    seq_along(1:nrow(data)), 
    function(i) {
        serialize(
            make_protobuf_feed(data[i, ]),
            connection = sprintf("%s/trip_update_%d.pb", dir, data$timestamp[i])
        )
    }
))

system("cp simulations/mock_server.js scripts/simulation_traveltime_01/")
system("ln -sf ../../simulations/node_modules scripts/simulation_traveltime_01/node_modules")
system("killall node")
system("cd scripts/simulation_traveltime_01 && nohup node mock_server.js > /dev/null 2>&1 &")

## --- run model on that data to estimate travel times 
library(transitr)
library(magrittr)
setwd(simdir)

if (dir.exists("etas")) unlink("etas", recursive = TRUE, force = TRUE)
dir.create("etas")

if (dir.exists("modeleval")) unlink("modeleval", recursive = TRUE, force = TRUE)
dir.create("modeleval")

if (file.exists("segment_states.csv")) unlink("segment_states.csv")
if (file.exists("segment_observations.csv")) unlink("segment_observations.csv")

if (dir.exists("history")) unlink("history", recursive = TRUE, force = TRUE)
dir.create("history")

# config
nw <- load_gtfs("../../at_gtfs.db", output = "etas.pb") %>%
    realtime_feed(sprintf("http://localhost:3000/sim_tt_01/trip_updates_only"),
                  response = "protobuf") %>%
    set_parameters(
        n_particles = 10000,
        system_noise = 1.5
    )

## reset server ID
RCurl::getURL(sprintf("localhost:3000/%s/reset", "sim_tt_01"))

model(nw)

## --- examine estimation accuracy
segfiles <- list.files("history", pattern = "segment_", full = T)
segtt <- lapply(segfiles, read_csv, 
    col_types = "iinn", 
    col_names = c("segment_id", "timestamp", "travel_time", "error")) %>%
    bind_rows() %>%
    mutate(truth = travel_times[-1])

# need nodes
con <- dbConnect(SQLite(), nw$database)
nodes <- con %>% tbl("shape_nodes") %>% filter(shape_id==!!shape$shape_id[1])
segs <- con %>% tbl("road_segments") %>% 
    filter(road_segment_id %in% !!segtt$segment_id)

segdata <- inner_join(nodes, segs, by = c("node_id" = "node_from")) %>%
    inner_join(segtt, by = c("road_segment_id" = "segment_id"), copy = TRUE) %>%
    arrange(node_sequence) %>% 
    select(road_segment_id, travel_time, error, truth)
