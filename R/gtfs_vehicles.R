create_vehicles <- function(db) {
    con <- db_connect(db)
    on.exit(db_close(con))
    if (RSQLite::dbExistsTable(con, "vehicles")) {
        stop("Vehicles table already exists")
    }

    res <- RSQLite::dbSendQuery(con,
               paste_nl(
                   "CREATE TABLE vehicles (",
                   "  vehicle_id TEXT PRIMARY KEY,",
                   "  trip_id TEXT,",
                   "  timestamp TEXT,",
                   "  position_latitude FLOAT,",
                   "  position_longitude FLOAT,",
                   "  distance FLOAT,",
                   "  speed FLOAT,",
                   "  progress INTEGER",
                   ")"))
    RSQLite::dbClearResult(res)
    res <- RSQLite::dbSendQuery(con, "PRAGMA journal_mode=wal")
    RSQLite::dbClearResult(res)
}

check_vehicles <- function(db) {
    con <- db_connect(db)
    on.exit(db_close(con))
    res <- identical(RSQLite::dbGetQuery(con, "PRAGMA table_info(vehicles)"),
              data.frame(cid = 0:7,
                         name = c("vehicle_id", "trip_id", "timestamp", "position_latitude",
                                  "position_longitude", "distance", "speed", "progress"),
                         type = c("TEXT", "TEXT", "TEXT", "FLOAT", "FLOAT", "FLOAT", "FLOAT", "INTEGER"),
                         notnull = 0L, dflt_value = as.logical(NA),
                         pk = c(1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                         stringsAsFactors = FALSE))
    res
}
