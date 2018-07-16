create_stops <- function(object) {
    con <- object$connection
    if (RSQLite::dbExistsTable(con, "stops")) {
        stop("Stops table already exists")
    }

    res <- dbSendQuery(
        con,
        paste_nl(
            "CREATE TABLE stops (",
            "  stop_id TEXT PRIMARY KEY,",
            "  stop_lat DOUBLE,",
            "  stop_lon DOUBLE,",
            "  stop_code TEXT,",
            "  stop_name TEXT,",
            "  stop_desc TEXT,",
            "  zone_id TEXT,",
            "  parent_station TEXT,",
            "  location_type INTEGER,",
            "  version DOUBLE",
            ")"))
    dbClearResult(res)
}

update_stops <- function(object, file) {
    existing <- dbGetQuery(object$connection,
                           "SELECT stop_id FROM stops")
    tbl <- read.csv(file, header = TRUE)
    stops <- data.frame(stop_id = as.character(tbl$stop_id),
                        stop_lat = as.numeric(tbl$stop_lat),
                        stop_lon = as.numeric(tbl$stop_lon),
                        stop_code = as.character(tbl$stop_code),
                        stop_name = as.character(tbl$stop_name),
                        stop_desc = as.character(tbl$stop_desc),
                        zone_id = as.character(tbl$zone_id),
                        parent_station = as.character(tbl$parent_station),
                        location_type = as.integer(tbl$location_type),
                        version = as.numeric(gsub(".+_v", "", tbl$stop_id)))
    stops <- stops[!stops$stop_id %in% existing, ]
    dbWriteTable(object$connection, "stops", stops, append = TRUE)
}
