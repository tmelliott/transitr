library(tidyverse)
library(RProtoBuf)
library(RSQLite)
library(dbplyr)
library(pbapply)
readProtoFiles("src/vendor/protobuf/gtfs-realtime.proto")

## Fetch the archive from the r-pi
if (FALSE) {
    date <- as.POSIXct("2018-08-24")
    system(sprintf("scp tom@130.216.51.230:/mnt/storage/history/%s/%s/%s/archive_%s.zip archive.zip",
                   format(date, "%Y"), format(date, "%m"), format(date, "%d"),
                   format(date, "%Y_%m_%d")))
}

## Extract the contents and read into memory
files <- unzip("archive.zip", exdir = tempdir())

vfiles <- files[grepl("vehicle_locations", files)]
vps <- do.call(bind_rows, pblapply(vfiles[1:100], function(f) {
    feed <- read(transit_realtime.FeedMessage, f)
    do.call(bind_rows, lapply(feed$entity, function(x) {
        if (!x$has('vehicle') ||
            !x$vehicle$has('vehicle') ||
            !x$vehicle$has('trip') ||
            !x$vehicle$has('position')) return(NULL)
        tibble(vehicle_id = x$vehicle$vehicle$id,
               timestamp = x$vehicle$timestamp,
               trip_id = x$vehicle$trip$trip_id,
               route_id = x$vehicle$trip$route_id,
               lon = x$vehicle$position$longitude,
               lat = x$vehicle$position$latitude)
    }))
}))

ggplot(vps, aes(lon, lat)) +
    geom_point() +
    coord_fixed(1.2)

## calculate distance to shape
vps$dist2line <- numeric(nrow(vps))
for (tid in unique(vps$trip_id)) {
    con <- dbConnect(SQLite(), "fulldata.db")
    sid <- con %>% tbl("trips") %>% filter(trip_id == tid) %>% select(shape_id) %>% collect %>% pluck("shape_id")
    shape <- con %>% tbl("shapes") %>% filter(shape_id == sid) %>% arrange(shape_pt_sequence) %>% collect
    dbDisconnect(con)
    
    sgeom <- sf::st_linestring(shape %>% select(shape_pt_lon, shape_pt_lat) %>% as.matrix)
    vps$dist2line[vps$trip_id == tid] <- 
        geosphere::dist2Line(vps[vps$trip_id == tid, c("lon", "lat")], sgeom)[,1]

    # p <- ggplot(shape, aes(shape_pt_lon, shape_pt_lat)) + geom_path() + coord_fixed(1.6) +
    #     geom_label(aes(lon, lat, label = round(dist2line)), data = vps[vps$trip_id == tid, ])
    #     # geom_point(aes(lon, lat), data = vps[vps$trip_id == tid, ])
    # print(p)
    # grid::grid.locator()
}

d <- vps$dist2line[vps$dist2line < 9 & vps$dist2line >= 0]

sig2 <- 3^2
hist(d^2 / sig2, 100, freq=F)
curve(dexp(x, 2), 0, 15, col="red", add=T)

xx <- rnorm(1000)^2 + rnorm(1000)^2
hist(xx, 100, freq=F)
curve(dexp(x, 0.5), 0, 15, col="red", add=T)