create_calendar <- function(con) {
    if (RSQLite::dbExistsTable(con, "calendar")) {
        stop("Calendar table already exists")
    }

    res <- RSQLite::dbSendQuery(
        con,
        paste_nl(
            "CREATE TABLE calendar (",
            "  service_id TEXT PRIMARY KEY,",
            "  monday INTEGER,",
            "  tuesday INTEGER,",
            "  wednesday INTEGER,",
            "  thursday INTEGER,",
            "  friday INTEGER,",
            "  saturday INTEGER,",
            "  sunday INTEGER,",
            "  start_date TEXT,",
            "  end_date TEXT,",
            "  version DOUBLE",
            ")"))
    RSQLite::dbClearResult(res)
}

update_calendar <- function(object, file) {
    existing <- RSQLite::dbGetQuery(object$connection,
                           "SELECT service_id FROM calendar")
    tbl <- read.csv(file, header = TRUE)
    calendar <- data.frame(service_id = as.character(tbl$service_id),
                           monday = as.integer(tbl$monday),
                           tuesday = as.integer(tbl$tuesday),
                           wednesday = as.integer(tbl$wednesday),
                           thursday = as.integer(tbl$thursday),
                           friday = as.integer(tbl$friday),
                           saturday = as.integer(tbl$saturday),
                           sunday = as.integer(tbl$sunday),
                           version = as.numeric(gsub(".+_v", "", tbl$service_id)))
    calendar <- calendar[!calendar$service_id %in% existing, ]
    RSQLite::dbWriteTable(object$connection, "calendar", calendar, append = TRUE)
}

check_calendar <- function(con) {
    identical(RSQLite::dbGetQuery(con, "PRAGMA table_info(calendar)"),
              data.frame(cid = 0:10,
                         name = c("service_id", "monday", "tuesday", "wednesday",
                                  "thursday", "friday", "saturday", "sunday",
                                  "start_date", "end_date", "version"),
                         type = c("TEXT", "INTEGER", "INTEGER", "INTEGER", "INTEGER",
                                  "INTEGER", "INTEGER", "INTEGER", "TEXT", "TEXT",
                                  "DOUBLE"),
                         notnull = 0L, dflt_value = as.logical(NA),
                         pk = c(1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                         stringsAsFactors = FALSE))
}
