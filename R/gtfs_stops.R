create_stops <- function(db) {
    con <- db_connect(db)
    if (RSQLite::dbExistsTable(con, "stops")) {
        stop("Stops table already exists")
    }

    res <- RSQLite::dbSendQuery(
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
    RSQLite::dbClearResult(res)

    db_close(con)
}

update_stops <- function(object, file) {
    con <- db_connect(object$database)
    existing <- RSQLite::dbGetQuery(con, "SELECT stop_id FROM stops")
    tbl <- utils::read.csv(file, header = TRUE)
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
    RSQLite::dbWriteTable(con, "stops", stops, append = TRUE)
    db_close(con)
}

check_stops <- function(db) {
    con <- db_connect(db)
    res <- identical(RSQLite::dbGetQuery(con, "PRAGMA table_info(stops)"),
              data.frame(cid = 0:9,
                         name = c("stop_id", "stop_lat", "stop_lon", "stop_code",
                                  "stop_name", "stop_desc", "zone_id",
                                  "parent_station", "location_type", "version"),
                         type = c("TEXT", "DOUBLE", "DOUBLE", "TEXT", "TEXT",
                                  "TEXT", "TEXT", "TEXT", "INTEGER", "DOUBLE"),
                         notnull = 0L, dflt_value = as.logical(NA),
                         pk = c(1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                         stringsAsFactors = FALSE))
    db_close(con)
    res
}
