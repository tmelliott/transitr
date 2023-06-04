create_calendar_dates <- function(db) {
    con <- db_connect(db)
    on.exit(db_close(con))
    if (RSQLite::dbExistsTable(con, "calendar_dates")) {
        stop("Calendar_Dates table already exists")
    }

    res <- RSQLite::dbSendQuery(
        con,
        paste_nl(
            "CREATE TABLE calendar_dates (",
            "  service_id TEXT,",
            "  date TEXT,",
            "  exception_type INTEGER,",
            "  PRIMARY KEY (service_id, date)",
            ")"))
    RSQLite::dbClearResult(res)
}

update_calendar_dates <- function(object, file) {
    con <- db_connect(object$database)
    on.exit(db_close(con))
    existing <- RSQLite::dbGetQuery(con, "SELECT service_id || date FROM calendar_dates")[[1]]
    tbl <- utils::read.csv(file, header = TRUE)
    calendar_dates <- data.frame(service_id = as.character(tbl$service_id),
                                 date = as.character(tbl$date),
                                 exception_type = as.integer(tbl$exception_type))
    calendar_dates <-
        calendar_dates[
            !paste0(calendar_dates$service_id, calendar_dates$date) %in% existing, ]
    RSQLite::dbWriteTable(con, "calendar_dates",
                          calendar_dates, append = TRUE)
}

check_calendar_dates <- function(db) {
    con <- db_connect(db)
    on.exit(db_close(con))
    res <- identical(RSQLite::dbGetQuery(con, "PRAGMA table_info(calendar_dates)"),
              data.frame(cid = 0:2,
                         name = c("service_id", "date", "exception_type"),
                         type = c("TEXT", "TEXT", "INTEGER"),
                         notnull = 0L, dflt_value = as.logical(NA),
                         pk = c(1L, 2L, 0L),
                         stringsAsFactors = FALSE))
    res
}
