create_versions <- function(con) {
    if (RSQLite::dbExistsTable(con, "versions")) {
        stop("Versions table already exists")
    }

    res <- RSQLite::dbSendQuery(
        con,
        paste_nl(
            "CREATE TABLE versions (",
            "  version_string TEXT PRIMARY KEY,",
            "  version DOUBLE,",
            "  startdate TEXT,",
            "  enddate TEXT",
            ")"))
    RSQLite::dbClearResult(res)
}

update_versions <- function(object) {
    if (!has_version_api(object)) return()
    
    r <- send(object$apis$version)
    if (r$status != 200) {
        warning("Unable to fetch Version API")
        return()
    }
    v <- jsonlite::fromJSON(httr::content(r, "text"))$response
    v$version_string <- v$version
    v$version <- as.numeric(gsub(".+_v", "", v$version))
    v$startdate <- format(as.POSIXct(v$startdate), "%Y%m%d")
    v$enddate <- format(as.POSIXct(v$enddate), "%Y%m%d")
    RSQLite::dbWriteTable(object$connection, "versions", v, append = TRUE)
    cat("Versions updated\n")
}

version_api <- function(object, url, headers = list()) {
    object$apis$version <- api(url, headers)
    object
}
