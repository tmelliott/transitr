##' Create a new GTFS database, optionally from a path.
##'
##' The GTFS data, if specified, should be specifided as a series of text files.
##' See \code{update} for more info.
##' @title Create GTFS database
##' @param from where the GTFS data comes from
##' @param to where the GTFS data is going to be stored
##' @return a \code{trgtfs} object
##' @author Tom Elliott
##' @export
create_gtfs <- function(from, to = tempfile()) {
    if (!requireNamespace("RSQLite", quietly = TRUE)) {
        stop("Please install the `RSQLite` package first,\n",
             "  or connect to a database manually and use `load_gtfs` instead.")
    }
    con <- RSQLite::dbConnect(RSQLite::SQLite(), to)
    create_tables(con)
    nw <- load_gtfs(con)
    if (!missing(from)) {
        update(nw, from)
    }
    nw
}

create_tables <- function(con) {
    lapply(gtfs_tables(), function(tbl) {
        eval(parse(text = sprintf("create_%s", tbl)))(con)
    })
}

check_tables <- function(con) {
    res <- sapply(gtfs_tables(), function(tbl) {
        eval(parse(text = sprintf("check_%s", tbl)))(con)
    })
    all(res)
}

##' Load an existing GTFS database into R for use with transitr.
##'
##' @title Load GTFS database
##' @param con a database connection from \code{dbConnect}
##' @return a \code{trgtfs} object
##' @author Tom Elliott
##' @export
load_gtfs <- function(con) {
    if (!check_tables(con)) {
        stop("Oops, some of the tables aren't right...")
    }
    structure(list(connection = con),
              class = "trgtfs")
}

##' Update a GTFS database object with new data.
##'
##' The data is expected to consist of a series of GTFS text/csv files,
##' and can be specified as either
##' \begin{itemize}
##'   \item a URL to zipped GTFS data
##'   \item a path to a zipped GTFS file
##'   \item a path to a GTFS directory of text files
##' \end{itemize}
##' @title Update GTFS data
##' @param object the \code{trgtfs} database object to update
##' @param from where the new data should come from
##' @return a \code{trgtfs} object (with updated data)
##' @author Tom Elliott
##' @export
update.trgtfs <- function(object, from) {
    if (dir.exists(from)) {
        cat("Create from directory\n")
        fn <- "dir"
    } else if (file.exists(from)) {
        cat("Create from ZIP\n")
        fn <- "zip"
    } else {
        cat("Create from URL\n")
        fn <- "url"
    }

    ## Create a function call e.g., `.update_url(object, from)`
    eval(parse(text = sprintf(".update_%s", fn)))(object, from)
}

### Update methods for various types of data location
.update_dir <- function(object, dir) {
    lapply(list.files(dir, full.names = TRUE), function(file) {
        table <- get_table_name(file)
        if (check_valid_table(table)) {
            update_table(object, table, file)
        } else {
            cat("Skipping", file, "\n")
        }
    })
}

.update_zip <- function(object, file) {
    d <- file.path(tempdir(), "data")
    unzip(file, exdir = d)
    .update_dir(object, d)
}

.update_url <- function(object, url) {
    f <- tempfile(fileext = ".zip")
    download.file(url, f)
    .update_zip(object, f)
}

### Update methods for different tables
update_table <- function(object, table, file) {
    eval(parse(text = sprintf("update_%s", table)))(object, file)
}

### Some helper functions

get_table_name <- function(file) {
    ext <- tools::file_ext(file)
    if (!ext %in% c('txt', 'csv')) {
        warning("Unrecognised extension: ",
                "transitr only knows how to handle txt and csv files")
    }
    tools::file_path_sans_ext(basename(file))
}

gtfs_tables <- function() {
    c("agency", "routes", "trips", "shapes", "stops",
      "stop_times", "calendar", "calendar_dates")
}

check_valid_table <- function(name) {
    name %in% gtfs_tables()
}
