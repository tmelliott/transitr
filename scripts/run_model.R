library(transitr)

library(magrittr)

url <- "https://cdn01.at.govt.nz/data/gtfs.zip"
outdated <- function(url) {
    date <- gsub("Last-Modified: |\r", "", system(sprintf("curl -sI %s | grep Last-Modified", url), intern = TRUE))
    last.mod <- as.POSIXct(date, tz = "GMT", format = "%a, %d %b %Y %H:%M:%S GMT")

    cur.mod <- ".gtfs_date"
    if (!file.exists(cur.mod)) {
        writeLines(as.character(last.mod), cur.mod)
        return(TRUE)
    }
    prev.mod <- as.POSIXct(readLines(cur.mod), tz = "GMT")
    if (prev.mod < last.mod) {
        writeLines(as.character(last.mod), cur.mod)
        return(TRUE)
    }
    FALSE
}

if (!file.exists("fulldata.db")) {
    nw <- create_gtfs(url, db = "fulldata.db", output = "at_predictions.pb") %>% construct()
} else {
    nw <- load_gtfs("fulldata.db", output = "at_predictions.pb")
    if (outdated(url)) {
        try(nw %>% update(url))
        nw %>% construct()
    }
}

nw <- nw %>%
    realtime_feed(
        c("https://dl.dropboxusercontent.com/s/1fvto9ex649mkri/vehicle_locations.pb?dl=1",
          "https://dl.dropboxusercontent.com/s/v680zazqjgm3whs/trip_updates.pb?dl=1"),
        response = "protobuf")

if (file.exists("config.json")) {
    config <- jsonlite::read_json("config.json")
    nw$output <- sprintf("at_predictions.pb")
    nw <- do.call(set_parameters, c(list(nw), config))
    nw %>% model
} else {
    ## run with defaults
    cat(" +++ No config.json file found - running with defaults +++\n")
    nw %>% set_parameters(n_core = 2) %>% model 
}

