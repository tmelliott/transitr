url <- "https://cdn01.at.govt.nz/data/gtfs.zip"
fzip <- tempfile(fileext = ".zip")
fdir <- file.path(tempdir(), "data")

download.file(url, fzip, quiet = TRUE)
unzip(fzip, exdir = fdir)

ROUTES <- c('881', '27H', '27W', '25L', '25B', '70', '75')
VERSION <- '67.28'

routes <- read.csv(file.path(fdir, 'routes.txt'))
routes <- routes[routes$route_short_name %in% ROUTES,]
routeids <- routes$route_id
routeids <- as.character(routeids[grepl(VERSION, routeids)])

trips <- read.csv(file.path(fdir, 'trips.txt'), stringsAsFactors = FALSE)
trips <- trips[trips$route_id %in% routeids,]
tripids <- tapply(trips$trip_id, trips$route_id, function(x) x[1])
trips <- trips[trips$trip_id %in% tripids, ]

shapes <- read.csv(file.path(fdir, 'shapes.txt'))
shapes <- shapes[shapes$shape_id %in% unique(trips$shape_id), ]

## with(shapes,
##      plot(shape_pt_lon, shape_pt_lat,
##           col = as.numeric(shape_id),
##           pch = 19, cex = 0.2, asp = 1.2))

stop_times <- read.csv(file.path(fdir, 'stop_times.txt'))
stop_times <- stop_times[stop_times$trip_id %in% tripids, ]
stop_times$stop_id <- gsub("-moved", "", stop_times$stop_id)

stops <- read.csv(file.path(fdir, 'stops.txt'))
stops <- stops[stops$stop_id %in% stop_times$stop_id, ]

calendar <- read.csv(file.path(fdir, 'calendar.txt'))
calendar <- calendar[calendar$service_id %in% trips$service_id, ]

calendar_dates <- read.csv(file.path(fdir, 'calendar_dates.txt'))
calendar_dates <- calendar_dates[1:10, ]

## rewrite everything
write.csv(routes, file.path(fdir, 'route.txt'), quote = FALSE)
write.csv(trips, file.path(fdir, 'trips.txt'), quote = FALSE)
write.csv(shapes, file.path(fdir, 'shapes.txt'), quote = FALSE)
write.csv(stops, file.path(fdir, 'stops.txt'), quote = FALSE)
write.csv(stop_times, file.path(fdir, 'stop_times.txt'), quote = FALSE)
write.csv(calendar, file.path(fdir, 'calendar.txt'), quote = FALSE)
write.csv(calendar_dates, file.path(fdir, 'calendar_dates.txt'), quote = FALSE)

d <- getwd()
setwd(fdir)
path <- file.path(system.file('inst', 'extdata', package = 'transitr'),
                  'auckland_gtfs.zip')
unlink(path)
zip(path, list.files())
setwd(d)

