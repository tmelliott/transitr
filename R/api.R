api <- function(url, headers) {
    structure(list(url = url, headers = headers),
              class = "trapi")
}

send <- function(api) {
    if (!inherits(api, "trapi"))
        stop("Not an API object")
    httr::GET(api$url, do.call(httr::add_headers, api$headers))
}
