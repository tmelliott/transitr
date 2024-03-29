create_stops <- function(db) {
    con <- db_connect(db)
    on.exit(db_close(con))
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
            "  node_id INTEGER,",
            "  version DOUBLE",
            ")"))
    RSQLite::dbClearResult(res)
}

update_stops <- function(object, file) {
    con <- db_connect(object$database)
    on.exit(db_close(con))
    existing <- RSQLite::dbGetQuery(con, "SELECT stop_id FROM stops")[[1]]
    tbl <- utils::read.csv(file, header = TRUE)
    stops <- data.frame(
        stop_id = as.character(tbl$stop_id),
        stop_lat = as.numeric(tbl$stop_lat),
        stop_lon = as.numeric(tbl$stop_lon),
        stop_code = as.character(tbl$stop_code),
        stop_name = as.character(tbl$stop_name),
        stop_desc = as.character(tbl$stop_desc),
        zone_id = as.character(tbl$zone_id),
        parent_station = as.character(tbl$parent_station),
        location_type = as.integer(tbl$location_type),
        node_id = NA,
        version = as.numeric(gsub(".+_v|-moved$", "", tbl$stop_id))
    )
    stops <- stops[!stops$stop_id %in% existing, ]
    RSQLite::dbWriteTable(con, "stops", stops, append = TRUE)
}

check_stops <- function(db) {
    con <- db_connect(db)
    on.exit(db_close(con))
    res <- identical(
        RSQLite::dbGetQuery(con, "PRAGMA table_info(stops)"),
        data.frame(
            cid = 0:10,
            name = c("stop_id", "stop_lat", "stop_lon", "stop_code",
                "stop_name", "stop_desc", "zone_id",
                "parent_station", "location_type", "node_id", "version"),
            type = c("TEXT", "DOUBLE", "DOUBLE", "TEXT", "TEXT",
                "TEXT", "TEXT", "TEXT", "INTEGER", "INTEGER", "DOUBLE"),
            notnull = 0L, 
            dflt_value = as.logical(NA),
            pk = c(1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
            stringsAsFactors = FALSE
        )
    )
    res
}
