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
