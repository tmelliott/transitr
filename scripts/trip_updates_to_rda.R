# convert trip_updates archive into an rda (containing a tibble with date and time separate)

exit <- function(x) {
    cat(x, "\n")
    if (!interactive()) quit("no")
}

if (!requireNamespace("ssh", quietly = TRUE)) exit("Please install the `ssh` package")

ca <- commandArgs(TRUE)
if (length(ca) == 0) exit("Please specify a date using `--args YYYY-MM-DD`")

DATE <- tryCatch({
        d <- as.Date(ca[1])
        d
    },
    error = function(e) exit("Invalid date. Please specify using `YYYY-MM-DD`"),
    finally = cat(" * processing data for", format(d, '%Y-%m-%d'), "\n")
)

PI <- Sys.getenv('pi_ip')
if (PI == "") exit("Please set the `pi_ip` environment variable")

cat(" * connecting to remote host", PI, "... ")
session <- ssh::ssh_connect(paste("tom", PI, sep = "@"))
cat("done\n")

## download the archive
dir <- file.path("/mnt", "storage", "history", format(DATE, "%Y/%m/%d"))
file <- sprintf("archive_%s.zip", format(DATE, "%Y_%m_%d"))


cat(" * downloading archive ...")
d <- tempdir()
arch <- ssh::scp_download(session, file.path(dir, file), to = d, verbose = FALSE)
cat("done\n")
f <- file.path(d, file)

files <- unzip(f, list = TRUE)$Name
tufiles <- files[grepl("trip_", files)]
cat(" * archive contains", length(tufiles), "trip updates\n")

cat(" * reading protobuf file ...")
library(RProtoBuf)
curd <- setwd("src/vendor/protobuf")
readProtoFiles("gtfs-realtime.proto")
setwd(curd)
cat(" done\n")

suppressPackageStartupMessages(library(tidyverse))
cat(" * reading trip updates ...")
rda <- sprintf("trip_updates_%s.rda", format(date, "%Y-%m-%d"))
pbapply::pblapply(
    tufiles, function(file) {
        pb <- unzip(f, file = file)
        feed <- read(transit_realtime.FeedMessage, pb)
        unlink(pb)
        if (length(feed$entity) == 0) return(NULL)
        lapply(feed$entity, function(e) {
            stus <- e$trip_update$stop_time_update
            if (length(stus) == 0) return(NULL)
            xdf <- tibble(
                vehicle_id = e$trip_update$vehicle$id,
                trip_id = e$trip_update$trip$trip_id,
                route_id = e$trip_update$trip$route_id,
                timestamp = as.POSIXct(e$trip_update$timestamp, origin = "1970-01-01"),
                stop_sequence = sapply(stus, function(stu) 
                    if (stu$has('stop_sequence')) stu$stop_sequence else NA
                ),
                type = sapply(stus, function(stu) ifelse(stu$has('arrival'), 'arrival', 'departure')),
                time = sapply(stus, function(stu) if (stu$has('arrival')) stu$arrival$time else stu$departure$time)
            )
        }) %>% bind_rows
    }
) %>% 
bind_rows %>%
save(file = rda)

cat("done\n * saved to", rda, "\n\n")
