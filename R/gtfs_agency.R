create_agency <- function(con) {
    if (RSQLite::dbExistsTable(con, "agency")) {
        stop("Agency table already exists")
    }

    res <- dbSendQuery(con,
               paste_nl(
                   "CREATE TABLE agency (",
                   "  agency_id TEXT PRIMARY KEY,",
                   "  agency_name TEXT,",
                   "  agency_url TEXT,",
                   "  agency_phone TEXT,",
                   "  agency_timezone TEXT,",
                   "  agency_lang TEXT",
                   ")"))
    dbClearResult(res)
}

update_agency <- function(object, file) {
    existing <- dbGetQuery(object$connection,
                           "SELECT agency_id FROM agency")
    tbl <- read.csv(file, header = TRUE)    
    agency <- data.frame(agency_id = as.character(tbl$agency_id),
                         agency_name = tbl$agency_name,
                         agency_url = as.character(tbl$agency_url),
                         agency_phone = as.character(tbl$agency_phone),
                         agency_lang = as.character(tbl$agency_lang),
                         agency_timezone = as.character(tbl$agency_timezone),
                         stringsAsFactors = FALSE)
    agency <- agency[!agency$agency_id %in% existing, ]
    dbWriteTable(object$connection, "agency", agency, append = TRUE)
}
