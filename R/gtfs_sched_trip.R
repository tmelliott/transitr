loadModule("TripModule", TRUE)

#' GTFS Trip constructor
#'
#' @param id Trip ID
#' @param route Route object, constructed using `route()`
#' @param direction_id Direction ID (0 or 1)
#' @param headsign Trip headsign
#' @param version GTFS Version
#' @return Rcpp_Trip object
#' @md
#' @export
trip <- function(id, route, direction_id = 0, headsign = "", version = 0) {
    if (!inherits(route, "Rcpp_Route")) {
        stop("route must be an Rcpp_Route object")
    }

    new("Rcpp_Trip", id, route, direction_id, headsign, version)
}

#' @describeIn trip Show method
#' @export
setMethod(
    "show", "Rcpp_Trip",
    function(object) {
        cat("<Trip>\n")
        cat("  trip_id:", object$trip_id(), "\n")
        cat(" route_id:", object$route()$route_id())
    }
)
