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
files <- list.files(file.path("simulations", "archive"), full.names = TRUE)

vfiles <- files[grepl("vehicle_locations", files)]
vps <- do.call(bind_rows, pblapply(vfiles, function(f) {
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




get_segment_data <- function(routes) {
    con <- dbConnect(SQLite(), "fulldata.db")
    on.exit(dbDisconnect(con))
    if (missing(routes)) {
        segments <- con %>% tbl("road_segments")
    } else {
        rids <- con %>% tbl("routes") %>%
            filter(route_short_name %in% routes &
                (route_long_name %like% "%To City%" |
                 route_long_name %like% "%To Britomart%")) %>%
            collect %>% pluck("route_id") %>% unique
        sids <- con %>% tbl("trips") %>%
            filter(route_id %in% rids) %>%
            collect %>% pluck("shape_id") %>% unique
        segids <- con %>% tbl("shape_segments") %>% 
            filter(shape_id %in% sids) %>%
            collect %>% pluck("road_segment_id") %>% unique
        segments <- con %>% tbl("road_segments") %>%
            filter(road_segment_id %in% segids)
    }
    intersections <- con %>% tbl("intersections")
    segments <- segments %>% 
        inner_join(intersections, by = c("int_from" = "intersection_id"), suffix = c("", "_start")) %>%
        inner_join(intersections, by = c("int_to" = "intersection_id"), suffix = c("", "_end")) %>%
        select(road_segment_id, length, intersection_lat, intersection_lon, intersection_lat_end, intersection_lon_end) %>%
        collect
    segments
}
sx <- get_segment_data() %>% filter(road_segment_id == "4691")
xr <- range(c(sx$intersection_lon, sx$intersection_lon_end))
yr <- range(c(sx$intersection_lat, sx$intersection_lat_end))
xr <- extendrange(xr, f = 1)
yr <- extendrange(yr, f = 1)

library(leaflet)

leaflet(vps %>% 
  filter(lon > xr[1] & lon < xr[2] & lat > yr[1] & lat < yr[2])
  ) %>%
  addTiles() %>%
  addCircles(~lon, ~lat, radius = 5, weight = 0)
