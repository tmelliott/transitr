create_stop_times <- function(db) {
    con <- db_connect(db)
    on.exit(db_close(con))
    if (RSQLite::dbExistsTable(con, "stop_times")) {
        stop("Stop_Times table already exists")
    }

    res <- RSQLite::dbSendQuery(
        con,
        paste_nl(
            "CREATE TABLE stop_times (",
            "  trip_id TEXT,",
            "  stop_id TEXT,",
            "  stop_sequence INTEGER,",
            "  arrival_time TEXT,",
            "  departure_time TEXT,",
            "  stop_headsign TEXT,",
            "  pickup_type INTEGER,",
            "  drop_off_type INTEGER,",
            "  shape_dist_traveled DOUBLE,",
            "  PRIMARY KEY (trip_id, stop_id, stop_sequence)",
            ")"))
    RSQLite::dbClearResult(res)
}

update_stop_times <- function(object, file) {
    con <- db_connect(object$database)
    on.exit(db_close(con))
    existing <- RSQLite::dbGetQuery(con, "SELECT trip_id || stop_id || stop_sequence FROM stop_times")[[1]]
    tbl <- utils::read.csv(file, header = TRUE)
    stop_times <- data.frame(trip_id = as.character(tbl$trip_id),
                             stop_id = as.character(tbl$stop_id),
                             stop_sequence = as.integer(tbl$stop_sequence),
                             arrival_time = as.character(tbl$arrival_time),
                             departure_time = as.character(tbl$departure_time),
                             stop_headsign = as.character(tbl$stop_headsign),
                             pickup_type = as.integer(tbl$pickup_type),
                             drop_off_type = as.integer(tbl$drop_off_type),
                             shape_dist_traveled = as.numeric(tbl$shape_dist_traveled)
                             )
    stop_times <-
        stop_times[
            !paste0(stop_times$trip_id, stop_times$stop_id,
                    stop_times$stop_sequence) %in% existing, ]
    RSQLite::dbWriteTable(con, "stop_times", stop_times, append = TRUE)
}

check_stop_times <- function(db) {
    con <- db_connect(db)
    on.exit(db_close(con))
    res <- identical(RSQLite::dbGetQuery(con, "PRAGMA table_info(stop_times)"),
              data.frame(cid = 0:8,
                         name = c("trip_id", "stop_id", "stop_sequence", "arrival_time",
                                  "departure_time", "stop_headsign", "pickup_type",
                                  "drop_off_type", "shape_dist_traveled"),
                         type = c("TEXT", "TEXT", "INTEGER", "TEXT", "TEXT", "TEXT",
                                  "INTEGER", "INTEGER", "DOUBLE"),
                         notnull = 0L, dflt_value = as.logical(NA),
                         pk = c(1L, 2L, 3L, 0L, 0L, 0L, 0L, 0L, 0L),
                         stringsAsFactors = FALSE))
    res
}
