##' Create a new GTFS database, optionally from a path.
##'
##' The GTFS data, if specified, should be specifided as a series of text files.
##' See \code{update} for more info.
##' @title Create GTFS database
##' @param source where the GTFS data comes from
##' @param db where the GTFS data is going to be stored
##' @param quiet logical, if \code{TRUE} progress output will be suppressed
##' @param output filename to output ETA predictions
##' @return a \code{trgtfs} object
##' @author Tom Elliott
##' @importFrom stats update
##' @export
create_gtfs <- function(source, db = tempfile(), quiet = FALSE, output = "predictions.pb") {
    if (!requireNamespace("RSQLite", quietly = TRUE)) {
        stop("Please install the `RSQLite` package first,\n",
             "  or connect to a database manually and use `load_gtfs` instead.")
    }
    
    create_tables(db)
    nw <- load_gtfs(db, output)
    if (!missing(source)) {
        update(nw, source, quiet = quiet)
    }
    nw
}

create_tables <- function(db) {
    lapply(gtfs_tables(), function(tbl) {
        eval(parse(text = sprintf("create_%s", tbl)))(db)
    })
    create_versions(db)
    create_vehicles(db)

    invisible(NULL)
}

check_tables <- function(db) {
    res <- sapply(c(gtfs_tables(), "vehicles"), function(tbl) {
        eval(parse(text = sprintf("check_%s", tbl)))(db)
    })
    all(res)
}

##' Load an existing GTFS database into R for use with transitr.
##'
##' @title Load GTFS database
##' @param db a database connection from \code{dbConnect}
##' @param output filename to output ETA predictions
##' @return a \code{trgtfs} object
##' @author Tom Elliott
##' @export
load_gtfs <- function(db, output = "predictions.pb") {
    if (!check_tables(db)) {
        stop("Oops, some of the tables aren't right...")
    }
    structure(list(database = db, apis = apis(), output = output,
                   parameters = list(n_core = 1L, 
                                     n_particles = 1000L, 
                                     gps_error = 5.0,
                                     system_noise = 1.0,
                                     pr_stop = 0.5,
                                     dwell_time = 10.0,
                                     gamma = 6.0,
                                     save_timings = FALSE)),
              class = "trgtfs")
}

##' Update a GTFS database object with new data.
##'
##' The data is expected to consist of a series of GTFS text/csv files,
##' and can be specified as either
##' \itemize{
##'   \item{}{a URL to zipped GTFS data}
##'   \item{}{a path to a zipped GTFS file}
##'   \item{}{a path to a GTFS directory of text files}
##' }
##' @title Update GTFS data
##' @param object the \code{trgtfs} database object to update
##' @param source where the new data should come from
##' @param quiet logical, if \code{TRUE} progress output will be suppressed
##' @param ... additional arguments, ignored
##' @return a \code{trgtfs} object (with updated data)
##' @author Tom Elliott
##' @export
update.trgtfs <- function(object, source, quiet = FALSE, ...) {
    if (dir.exists(source)) {
        fn <- "dir"
    } else if (file.exists(source)) {
        fn <- "zip"
    } else {
        fn <- "url"
    }

    ## Create a function call e.g., `.update_url(object, source)`
    eval(parse(text = sprintf(".update_%s", fn)))(object, source, quiet)
    ## Update versions
    update_versions(object)
}

### Update methods for various types of data location
.update_dir <- function(object, dir, quiet) {
    lapply(list.files(dir, full.names = TRUE), function(file) {
        table <- get_table_name(file)
        if (check_valid_table(table)) {
            if (!quiet) cat(" * creating", file, "\n")
            update_table(object, table, file)
        } else {
            if (!quiet) cat(" * skipping", file, "\n")
        }
    })
}

.update_zip <- function(object, file, quiet) {
    d <- file.path(tempdir(), "data")
    utils::unzip(file, exdir = d)
    .update_dir(object, d, quiet)
}

.update_url <- function(object, url, quiet) {
    f <- tempfile(fileext = ".zip")
    utils::download.file(url, f, quiet = quiet)
    .update_zip(object, f, quiet)
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

has_version_api <- function(object) {
    !is.null(object$apis$version)
}
