loadModule("VehicleModule", TRUE)

#' @export
vehicle <- function(id) {
    new("Rcpp_Vehicle", id)
}

setMethod(
    "show", "Rcpp_Vehicle",
    function(object) {
        id <- object$get_id()
        cat("Vehicle: ", id)
    }
)
