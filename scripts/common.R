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
                xdf <- tibble(
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
                qq <- stus[[length(stus)]]$getExtension(transit_network.eta)$quantiles
                if (length(qq) > 0) {
                    ## fetch quantiles and bind to xdf
                    quantiles <- sapply(qq, function(x) x$quantile)
                    qs <- sapply(stus, function(stu) {
                        qs <- sapply(stu$getExtension(transit_network.eta)$quantiles, 
                            function(x) ifelse(x$value == 0, NA, x$value))
                        if (length(qs) != length(quantiles)) return(integer(length(quantiles)))
                        qs
                    }) %>% t %>% as.tibble
                    names(qs) <- paste0("q", quantiles)
                    xdf <- bind_cols(xdf, qs)
                }
                xdf
            })
        ) %>%
            mutate(time = as.POSIXct(time, origin = "1970-01-01")) %>%
            group_by(trip_id)
        save(etas, file = rda)
    }
    etas
}

all_sims <- function(sim) {
    times <- gsub("etas_|\\.pb", "", list.files(file.path("simulations", sim, "etas"), ".pb")) %>%
        as.integer
    do.call(bind_rows, pbapply::pblapply(times, function(t) loadsim(sim, t)))
}

library(RSQLite)
library(dbplyr)
get_schedule <- function(trip) {
    con <- dbConnect(SQLite(), "fulldata.db")
    on.exit(dbDisconnect(con))
    con %>% tbl("stop_times") %>% filter(trip_id == trip) %>% arrange(stop_sequence) %>%
        select(trip_id, stop_sequence, arrival_time, departure_time) %>% collect
}
