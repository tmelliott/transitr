create_routes <- function(object) {
    con <- object$connection
    if (RSQLite::dbExistsTable(con, "routes")) {
        stop("Routes table already exists")
    }

    res <- dbSendQuery(
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
    dbClearResult(res)
}

update_routes <- function(object, file) {
    existing <- dbGetQuery(object$connection,
                           "SELECT route_id FROM routes")
    tbl <- read.csv(file, header = TRUE)
    routes <- data.frame(route_id = as.character(tbl$route_id),
                         route_short_name = as.character(tbl$route_short_name),
                         route_long_name = as.character(tbl$route_long_name),
                         route_type = as.integer(tbl$route_type),
                         agency_id = as.character(tbl$agency_id),
                         version = as.numeric(gsub(".+_v", "", tbl$route_id)))
    routes <- routes[!routes$route_id %in% existing, ]
    dbWriteTable(object$connection, "routes", routes, append = TRUE)
}
