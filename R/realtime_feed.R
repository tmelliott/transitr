
realtime_feed <- function(nw, url, headers = NULL, proto = FALSE) {
    if (proto) headers <- c(headers, with_headers('Accept' = 'application/x-protobuf'))
    nw$apis$realtime <- api(url, headers)
    nw
}