create_trips <- function(object) {
    con <- object$connection
    if (RSQLite::dbExistsTable(con, "trips")) {
        stop("Trips table already exists")
    }

    res <- dbSendQuery(
        con,
        paste_nl(
            "CREATE TABLE trips (",
            "  trip_id TEXT PRIMARY KEY,",
            "  route_id TEXT,",
            "  shape_id TEXT,",
            "  service_id TEXT,",
            "  block_id TEXT,",
            "  direction_id INTEGER,",
            "  trip_headsign TEXT,",
            "  version DOUBLE",
            ")"))
    dbClearResult(res)
}

update_trips <- function(object, file) {
    existing <- dbGetQuery(object$connection,
                           "SELECT trip_id FROM trips")
    tbl <- read.csv(file, header = TRUE)
    trips <- data.frame(trip_id = as.character(tbl$trip_id),
                        route_id = as.character(tbl$route_id),
                        shape_id = as.character(tbl$shape_id),
                        service_id = as.character(tbl$service_id),
                        block_id = as.character(tbl$block_id),
                        direction_id = as.integer(tbl$direction_id),
                        trip_headsign = as.character(tbl$trip_headsign),
                        version = as.numeric(gsub(".+_v", "", tbl$route_id)))
    trips <- trips[!trips$trip_id %in% existing, ]
    dbWriteTable(object$connection, "trips", trips, append = TRUE)
}
