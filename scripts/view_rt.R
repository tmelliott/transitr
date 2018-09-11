library(RSQLite)
library(ggplot2)
library(ggmap)
db <- "fulldata.db"

con <- dbConnect(SQLite(), db)
v <- dbGetQuery(con, 
    "SELECT timestamp, position_latitude as lat, position_longitude as lng, distance, speed*60*60/1000 as speed, progress FROM vehicles")
dbDisconnect(con)

v$timestamp <- as.POSIXct(as.integer(v$timestamp), origin = "1970-01-01")
v <- v[Sys.time() - v$timestamp < 120, ]

box <- c(174.4377, -37.32487, 175.0781, -36.54473)
map <- get_stamenmap(box, zoom = 11, maptype = "toner-lite")

p <- ggmap(map) + 
    geom_point(aes(lng, lat, colour = speed), data = v) +
    scale_colour_viridis_c(option="C")

dev.hold()
print(p)
dev.flush()







## something else
library(tidyverse)
library(RSQLite)
library(dbplyr)
library(ggmap)

doaplot <- function(id, which = 1, nrow = NULL, ncol = NULL) {
    v <- read.csv(sprintf("history/vehicle_%s.csv", vs[vn]), header = FALSE)
    p <- read.csv(sprintf("history/vehicle_%s_particles.csv", vs[vn]), header = FALSE)
    colnames(v) <- c("timestamp", "trip_id", "obs_lat", "obs_lon", "distance", "speed", "latitude", "longitude")
    colnames(p) <- c("timestamp", paste0("stop", 1:(ncol(p) - 1)))
    v <- v %>% as.tibble %>%
        filter(timestamp > min(timestamp) & trip_id == names(table(v$trip_id))[1]) %>%
        mutate(timestamp = as.POSIXct(timestamp, origin = "1970-01-1"))
    p <- p %>% as.tibble %>% 
        filter(timestamp %in% unique(v$timestamp)) %>%
        gather("stop_number", "arrival_time", -timestamp) %>% 
        filter(arrival_time > 0) %>% 
        mutate(stop_number = as.numeric(gsub("stop", "", stop_number)),
               timestamp = as.POSIXct(timestamp, origin = "1970-01-01"),
               arrival_time = pmin(Sys.time() + 3600, as.POSIXct(arrival_time, origin = "1970-01-01")))

    con <- dbConnect(SQLite(), "fulldata.db")
    tid <- v$trip_id[1]
    shape <- con %>% tbl("trips") %>% 
        filter(trip_id == tid) %>%
        inner_join(con %>% tbl("shapes"), "shape_id") %>%
        select(shape_pt_lat, shape_pt_lon) %>%
        arrange(shape_pt_sequence) %>%
        collect
    dbDisconnect(con)

    switch(which,
        {
            xr <- extendrange(shape$shape_pt_lon)
            yr <- extendrange(shape$shape_pt_lat)
            map <- get_stamenmap(c(xr[1], yr[1], xr[2], yr[2]), zoom = 14, maptype = "toner-lite")
            map %>% ggmap() + geom_path(aes(shape_pt_lon, shape_pt_lat), data = shape, lwd = 2, color = "steelblue") +
                geom_point(aes(x = longitude, y = latitude), data = v, color = "orange") + 
                geom_point(aes(obs_lon, obs_lat), data = v, color = "red", shape = 4) +
                facet_wrap(~timestamp, nrow = nrow, ncol = ncol)
        },
        p %>% ggplot(aes(timestamp, arrival_time)) + geom_point() + facet_wrap(~stop_number),
        p %>% ggplot(aes(arrival_time, stop_number)) + geom_point() + facet_wrap(~timestamp),
        p %>% ggplot(aes(arrival_time, stop_number)) + geom_point(aes(color = timestamp)),
        {
        }
    )
}

ps <- list.files("history", pattern = "*_particles.csv")
vs <- gsub("vehicle_|_particles.csv", "", ps)
vn <- 25
doaplot(vs[vn], 1, nrow = 4)
doaplot(vs[vn], 2)
doaplot(vs[vn], 3)
doaplot(vs[vn], 4)
