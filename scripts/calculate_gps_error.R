library(tidyverse)
library(RSQLite)
library(dbplyr)

con <- dbConnect(SQLite(), "fulldata.db")
vps <- con %>% tbl("vehicles") %>% filter(progress > 0 & progress < 100) %>% collect
dbDisconnect(con)

ggplot(vps, aes(position_longitude, position_latitude)) +
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
        geosphere::dist2Line(vps[vps$trip_id == tid, c("position_longitude", "position_latitude")], sgeom)[,1]

    p <- ggplot(shape, aes(shape_pt_lon, shape_pt_lat)) + geom_path() + coord_fixed(1.6) +
        geom_label(aes(position_longitude, position_latitude, label = round(dist2line)), data = vps[vps$trip_id == tid, ])
        # geom_point(aes(position_longitude, position_latitude), data = vps[vps$trip_id == tid, ])
    print(p)
    grid::grid.locator()
}

d <- vps$dist2line[vps$dist2line < 100 & vps$dist2line > 0]
hist(d, 100)
