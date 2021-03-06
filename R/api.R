#' Create an API call
#'
#' @param url the url to request
#' @param headers optional headers for the call
#' @param response the type of response expected
#' @return a \code{trapi} object
#' @author Tom Elliott
#' @export
api <- function(url, headers = NULL, response = "text") {
    structure(list(url = url, headers = headers, response = response),
              class = "trapi")
}

apis <- function() {
    structure(list(), class = "trapi.list")
}

#' @describeIn api Api Headers
#' @param ... headers converted into a list
#' @export
with_headers <- function(...) list(...)

send <- function(api) {
    if (!inherits(api, "trapi"))
        stop("Not an API object")
    httr::GET(api$url, do.call(httr::add_headers, api$headers))
}

print.trapi <- function(x, ...) {
    cat(x$url, "\n")
    if (!is.null(x$headers)) {
        for (i in seq_along(x$headers)) {
            name <- names(x$headers)[i]
            value <- x$headers[[i]]
            cat("  ", glue::glue("'{name}': '{value}'"), "\n")
        }
    }
}

print.trapi.list <- function(x, ...) {
    if (length(x) == 0) {
        cat("No APIs connected\n")
        return()
    }
    for (i in seq_along(x)) {
        cat(sep="", " ", glue::glue("{names(x)[i]} ({x[[i]]$response}): "))
        print(x[[i]])
    }
    cat("\n")
}