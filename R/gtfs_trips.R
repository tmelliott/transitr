create_trips <- function(con) {
    if (RSQLite::dbExistsTable(con, "trips")) {
        stop("Trips table already exists")
    }

    res <- RSQLite::dbSendQuery(
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
    RSQLite::dbClearResult(res)
}

update_trips <- function(object, file) {
    existing <- RSQLite::dbGetQuery(object$connection,
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
    RSQLite::dbWriteTable(object$connection, "trips", trips, append = TRUE)
}

check_trips <- function(con) {
    identical(RSQLite::dbGetQuery(con, "PRAGMA table_info(trips)"),
              data.frame(cid = 0:7,
                         name = c("trip_id", "route_id", "shape_id", "service_id",
                                  "block_id", "direction_id", "trip_headsign",
                                  "version"),
                         type = c("TEXT", "TEXT", "TEXT", "TEXT", "TEXT",
                                  "INTEGER", "TEXT", "DOUBLE"),
                         notnull = 0L, dflt_value = as.logical(NA),
                         pk = c(1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                         stringsAsFactors = FALSE))
}
