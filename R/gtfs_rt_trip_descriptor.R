loadModule("TripDescriptorModule", TRUE)

#' @export
trip_descriptor <- function(trip, direction_id = 0L, start_time = "",
                            start_date = "", schedule_relationship = "") {
    if (!is(trip, "Rcpp_Trip")) {
        stop("trip must be an Rcpp_Trip object")
    }

    new(
        "Rcpp_TripDescriptor",
        trip, direction_id, start_time, start_date, schedule_relationship
    )
}

setMethod(
    "show", "Rcpp_TripDescriptor",
    function(object) {
        cat("<TripDescriptor>\n")
        cat("  trip_id: ", object$trip()$trip_id(), "\n")
        cat("  route_id: ", object$route()$route_id())
    }
)
