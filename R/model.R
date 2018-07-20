model <- function(nw, n.particles = 500, cores = 1L) {
    if (!inherits(nw, "trgtfs")) stop("Not a GTFS transit object. See ?create_gtfs")

    if (!check_tables(nw)) stop("GTFS tables don't appear to be valid.")

    if (is.null(nw$apis$realtime)) stop("Please attach a realtime API feed, `realtime_feed()`")
    # x <- send(nw$apis$realtime)
    # if (x$status != 200) {
    #     print(x)
    #     stop("Unable to access the API feed - have you added the key?")
    # }

    run_realtime_model(nw, n.particles, cores)
}