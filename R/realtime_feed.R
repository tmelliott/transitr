##' Connect to a realtime GTFS feed
##'
##' @title Connect to GTFS feed
##' @param nw a \code{trgtfs} object
##' @param url the API url
##' @param headers a header object from \code{with_headers} defining any GET request headers (e.g., api keys)
##' @param response the type of response expected (json or protobuf)
##' @return a \code{trgtfs} object with feed info
##' @author Tom Elliott
##' @export
realtime_feed <- function(nw, url, headers = NULL, response = c('json', 'protobuf')) {
    response <- match.arg(response)
    if (response == "protobuf") 
        headers <- c(headers, with_headers('Accept' = 'application/x-protobuf'))
    nw$apis$realtime <- api(url, headers, response)
    nw
}
