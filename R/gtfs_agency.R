create_agency <- function(db) {
    con <- db_connect(db)
    if (RSQLite::dbExistsTable(con, "agency")) {
        stop("Agency table already exists")
    }

    res <- RSQLite::dbSendQuery(con,
               paste_nl(
                   "CREATE TABLE agency (",
                   "  agency_id TEXT PRIMARY KEY,",
                   "  agency_name TEXT,",
                   "  agency_url TEXT,",
                   "  agency_phone TEXT,",
                   "  agency_timezone TEXT,",
                   "  agency_lang TEXT",
                   ")"))
    RSQLite::dbClearResult(res)

    db_close(con)
}

update_agency <- function(object, file) {
    con <- db_connect(object$database)
    existing <- RSQLite::dbGetQuery(con, "SELECT agency_id FROM agency")
    tbl <- utils::read.csv(file, header = TRUE)    
    agency <- data.frame(agency_id = as.character(tbl$agency_id),
                         agency_name = tbl$agency_name,
                         agency_url = as.character(tbl$agency_url),
                         agency_phone = as.character(tbl$agency_phone),
                         agency_lang = as.character(tbl$agency_lang),
                         agency_timezone = as.character(tbl$agency_timezone),
                         stringsAsFactors = FALSE)
    agency <- agency[!agency$agency_id %in% existing, ]
    RSQLite::dbWriteTable(con, "agency", agency, append = TRUE)
    db_close(con)
}

check_agency <- function(db) {
    con <- db_connect(db)
    res <- identical(RSQLite::dbGetQuery(con, "PRAGMA table_info(agency)"),
              data.frame(cid = 0:5,
                         name = c("agency_id", "agency_name", "agency_url",
                                  "agency_phone", "agency_timezone", "agency_lang"),
                         type = "TEXT", notnull = 0L, dflt_value = as.logical(NA),
                         pk = c(1L, 0L, 0L, 0L, 0L, 0L),
                         stringsAsFactors = FALSE))
    db_close(con)
    res
}
