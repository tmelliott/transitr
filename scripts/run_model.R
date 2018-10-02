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
        nw %>% update(url) %>% construct()
    }
}

nw <- nw %>%
    realtime_feed("https://dl.dropboxusercontent.com/s/1fvto9ex649mkri/vehicle_locations.pb?dl=1", 
                  response = "protobuf")

if (Sys.info()["nodename"] == "certellprd01") {
    nw <- nw %>% set_parameters(n_particles = 5000, n_core = 6, gps_error = 5, system_noise = 0.5)
    if (Sys.getenv("GPSERROR") == "")
        nw %>% model
    else {
        nw$output <- sprintf("at_predictions_%s.pb", Sys.getenv("GPSERROR"))
        nw %>% set_parameters(gps_error = as.numeric(Sys.getenv("GPSERROR"))) %>% model
    }

} else {
    nw %>% set_parameters(n_particles = 1000, n_core = 1, gps_error = 5, system_noise = 0.5) %>% model
}


