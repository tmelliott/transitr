create_stop_times <- function(object) {
    con <- object$connection
    if (RSQLite::dbExistsTable(con, "stop_times")) {
        stop("Stop_Times table already exists")
    }

    res <- dbSendQuery(
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
    dbClearResult(res)
}

update_stop_times <- function(object, file) {
    existing <- dbGetQuery(object$connection,
                           "SELECT trip_id || stop_id || stop_sequence FROM stop_times")
    tbl <- read.csv(file, header = TRUE)
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
    dbWriteTable(object$connection, "stop_times", stop_times, append = TRUE)
}
