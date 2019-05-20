source("scripts/common.R")

## --- load a route and its stops (for a trip)
db <- "fulldata.db"
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
dbDisconnect(con)


## --- generate some simulated data with known travel times
# what we want ...
pb <- list.files("simulations/archive", full = TRUE)[1]
f <- read(transit_realtime.FeedMessage, pb)

Date <- "2019-04-01"
start <- as.integer(as.POSIXct(
    sprintf("%s %s", Date, stops$arrival_time[1])
))  
data <- data.frame(
    vehicle_id = "TEST",
    trip_id = trip$trip_id,
    route_id = route$route_id,
    timestamp = start,
    type = "arrival",
    stop_sequence = 1,
    stop_id = stops$stop_id[1],
    stringsAsFactors = FALSE
)
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
                            departure = new(transit_realtime.TripUpdate.StopTimeEvent,
                                time = data$time
                            )
                        )
                    )
                )
            )
        )
    )
}


## --- run model on that data to estimate travel times 



## --- examine estimation accuracy
