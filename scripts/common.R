library(tidyverse)
library(RProtoBuf)

curd <- setwd("src/vendor/protobuf")
readProtoFiles("gtfs-realtime-ext.proto")
setwd(curd)

loadsim <- function(sim, time) {
    time <- as.integer(time)
    pb <- file.path("simulations", sim, "etas", sprintf("etas_%s.pb", time))
    rda <- file.path("simulations", sim, "etas", sprintf("etas_%s.rda", time))
    if (!file.exists(pb)) {
        warning("That file doesn't exist:", pb)
        return(NULL)
    }
    if (file.exists(rda)) {
        load(rda)
    } else {
        feed <- read(transit_realtime.FeedMessage, pb)
        if (length(feed$entity) == 0) return(NULL)
        ent <- feed$entity
        etas <- do.call(bind_rows, 
            lapply(ent, function(e) {
                stus <- e$trip_update$stop_time_update
                if (length(stus) == 0) return(NULL)
                tibble(
                    vehicle_id = e$trip_update$vehicle$id,
                    trip_id = e$trip_update$trip$trip_id,
                    route_id = e$trip_update$trip$route_id,
                    timestamp = as.POSIXct(time, origin = "1970-01-01"),
                    stop_sequence = sapply(stus, function(stu) 
                        if (stu$has('stop_sequence')) stu$stop_sequence else NA
                        ),
                    time = sapply(stus, function(stu)
                        if (stu$getExtension(transit_network.eta)$has('estimate') &&
                            stu$getExtension(transit_network.eta)$estimate > 0) 
                            stu$getExtension(transit_network.eta)$estimate 
                        else NA
                        )
                )
            })
        ) %>%
            mutate(time = as.POSIXct(time, origin = "1970-01-01")) %>%
            group_by(trip_id)
        save(etas, file = rda)
    }
    etas
}

