loadModule("RouteModule", TRUE)

#' Gtfs 'Route' Constructor
#'
#' @param id Route ID
#' @param short_name Route short name
#' @param long_name Route long name
#' @param type Route type
#' @param agency_id Agency ID
#' @param version GTFS Version
#' @return Rcpp_Route object
#' @md
#' @export
route <- function(id, short_name = "", long_name = "", type = 0L,
                  agency_id = "", version = 0) {
    new(
        "Rcpp_Route",
        id, short_name, long_name, type, agency_id, version
    )
}

#' @include GenericMethods.R
#' @describeIn route Show method
#' @export
setMethod(
    "show", "Rcpp_Route",
    function(object) {
        cat("<Route>\n")
        cat("  route_id:", object$route_id())
    }
)
