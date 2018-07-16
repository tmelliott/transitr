create_routes <- function(con) {
    if (RSQLite::dbExistsTable(con, "routes")) {
        stop("Routes table already exists")
    }

    res <- RSQLite::dbSendQuery(
        con,
        paste_nl(
            "CREATE TABLE routes (",
            "  route_id TEXT PRIMARY KEY,",
            "  route_short_name TEXT,",
            "  route_long_name TEXT,",
            "  route_type INTEGER,",
            "  agency_id TEXT,",
            "  version DOUBLE",
            ")"))
    RSQLite::dbClearResult(res)
}

update_routes <- function(object, file) {
    existing <- RSQLite::dbGetQuery(object$connection,
                           "SELECT route_id FROM routes")
    tbl <- utils::read.csv(file, header = TRUE)
    routes <- data.frame(route_id = as.character(tbl$route_id),
                         route_short_name = as.character(tbl$route_short_name),
                         route_long_name = as.character(tbl$route_long_name),
                         route_type = as.integer(tbl$route_type),
                         agency_id = as.character(tbl$agency_id),
                         version = as.numeric(gsub(".+_v", "", tbl$route_id)))
    routes <- routes[!routes$route_id %in% existing, ]
    RSQLite::dbWriteTable(object$connection, "routes", routes, append = TRUE)
}

check_routes <- function(con) {
    identical(RSQLite::dbGetQuery(con, "PRAGMA table_info(routes)"),
              data.frame(cid = 0:5,
                         name = c("route_id", "route_short_name", "route_long_name",
                                  "route_type", "agency_id", "version"),
                         type = c("TEXT", "TEXT", "TEXT", "INTEGER", "TEXT", "DOUBLE"),
                         notnull = 0L, dflt_value = as.logical(NA),
                         pk = c(1L, 0L, 0L, 0L, 0L, 0L),
                         stringsAsFactors = FALSE))
}
