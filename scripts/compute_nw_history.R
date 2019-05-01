library(tidyverse)
library(RSQLite)
library(dbplyr)
library(RProtoBuf)

curd <- setwd("src/vendor/protobuf")
readProtoFiles("gtfs-realtime-ext.proto")
setwd(curd)

# 0. set up database
db <- "history.sqlite"
con <- dbConnect(SQLite(), db)

if (!dbExistsTable(con, "trip_updates")) {
    r <- dbSendQuery(con,
        paste(sep = "\n",
            "CREATE TABLE trip_updates (",
            "  date TEXT,",
            "  trip_id TEXT,",
            "  stop_id TEXT,",
            "  arrival_time INTEGER,",
            "  departure_time INTEGER,",
            "  dwell_time INTEGER",
            ")"
        )
    )
    dbClearResult(r)
    rm("r")
}

if (!dbExistsTable(con, "segments")) {
    r <- dbSendQuery(con,
        paste(sep = "\n",
            "CREATE TABLE segments (",
            "  date TEXT,",
            "  segment_id TEXT,",
            "  trip_id TEXT,",
            "  start_time INTEGER,",
            "  end_time INTEGER,",
            "  travel_time INTEGER",
            ")"
        )
    )
    dbClearResult(r)
    rm("r")
}

if ( ! db_has_table(con, "segments") )
    stop("Error connecting to database or creating segments table")

if ( ! dbConnect(SQLite(), db) %>% db_has_table("trip_updates") )
    stop("Error connecting to database or creating trip_updates table")

# 1. take a vector of dates
dates <- as.Date("2019-04-01")

for (d in seq_along(dates)) {
    date <- as.Date(dates[d])
    # 2. read trip updates
    archive_file <- sprintf("/mnt/storage/history/%s/archive_%s.zip",
        format(date, "%Y/%m/%d"),
        format(date, "%Y_%m_%d")
    )

    tmp <- tempfile()
    fetch <- system(
        sprintf("scp tom@%s:%s %s",
            Sys.getenv("pi_ip"),
            archive_file,
            tmp
        )
    )

    if (fetch != 0)
        stop("There was an error fetching the ZIP archive from the Pi.")

    # queries
    check_qry <- "SELECT count(*) as n FROM trip_updates WHERE date=? AND trip_id=? AND stop_id=?"
    insert_qry <- "INSERT INTO trip_updates VALUES (?, ?, ?, ?, ?, 0)"
    update_arr_qry <- "UPDATE trip_updates SET arrival_time=? WHERE date=? AND trip_id=? AND stop_id=?"
    update_dep_qry <- "UPDATE trip_updates SET departure_time=? WHERE date=? AND trip_id=? AND stop_id=?"

    all_files <- unzip(tmp, exdir = tempdir())
    tu_files <- all_files[grepl("trip_update", all_files)]
    con <- dbConnect(SQLite(), db)
    for (tu_file in tu_files) {
        tu_feed <- read(transit_realtime.FeedMessage, tu_file)
        for (i in seq_along(1:length(tu_feed$entity))) {
            e <- tu_feed$entity[[i]]
            # insert or update?
            q <- dbSendQuery(con, check_qry)
            dbBind(q,
                list(
                    date,
                    e$trip_update$trip$trip_id,
                    e$trip_update$stop_time_update[[1]]$stop_id
                )
            )
            res <- dbFetch(q)$n
            dbClearResult(q)
            if (res == 0) {
                # insert
                q <- dbSendQuery(con, insert_qry)
                dbBind(q,
                    list(
                        date,
                        e$trip_update$trip$trip_id,
                        e$trip_update$stop_time_update[[1]]$stop_id,
                        e$trip_update$stop_time_update[[1]]$arrival$time,
                        e$trip_update$stop_time_update[[1]]$departure$time
                    )
                )
                dbClearResult(q)
            } else {
                # update
                if (e$trip_update$stop_time_update[[1]]$arrival$time > 0) {
                    q <- dbSendQuery(con, update_arr_qry)
                    dbBind(q,
                        list(
                            e$trip_update$stop_time_update[[1]]$arrival$time,
                            date,
                            e$trip_update$trip$trip_id,
                            e$trip_update$stop_time_update[[1]]$stop_id
                        )
                    )
                    dbClearResult(q)
                } else {
                    q <- dbSendQuery(con, update_dep_qry)
                    dbBind(q,
                        list(
                            e$trip_update$stop_time_update[[1]]$departure$time,
                            date,
                            e$trip_update$trip$trip_id,
                            e$trip_update$stop_time_update[[1]]$stop_id
                        )
                    )
                    dbClearResult(q)
                }
            }
        }
    }

    unlink()

    # and compute dwell times in DB
    r <- dbSendQuery(con,
        "UPDATE trip_updates SET dwell_time = departure_time - arrival_time WHERE arrival_time > 0 AND departure_time > 0"
    )
    dbClearResult(r)
    # and clean up (negative dwell time, or 5+ min)
    r <- dbSendQuery(con,
        "DELETE FROM trip_updates WHERE dwell_time <= 0 OR dwell_time > 600"
    )
    dbClearResult(r)

    dbDisconnect(con)

    # 3. assign travel times to road segments
    # 3a. get trips that exist
    q <- dbSendQuery(con,
        "SELECT DISTINCT trip_id FROM trip_updates WHERE date=?"
    )
    dbBind(q, list(date))
    trips <- dbFetch(q)$trip_id
    dbClearResult(q)

    # 3b. create segments for stop sequences
    for (trip in trips) {
        con2 <- dbConnect(SQLite(), "fulldata.db")
        q <- dbSendQuery(con2,
            "SELECT stop_id FROM stop_times WHERE trip_id=? ORDER BY stop_sequence"
        )
        dbBind(q, list(trip))
        stops <- dbFetch(q)$stop_id
        dbClearResult(q)
        dbDisconnect(con2)

        sids <- gsub("-.*", "", stops)
        segment_ids <- paste(sids[-length(sids)], sids[-1], sep = ":")
        for (i in seq_along(segment_ids)) {
            q <- dbSendQuery(con,
                "SELECT count(*) as n FROM segments WHERE date=? AND segment_id=? AND trip_id=?"
            )
            dbBind(q,
                list(
                    date,
                    segment_ids[i],
                    trip
                )
            )
            res <- dbFetch(q)$n
            dbClearResult(q)
            if (res > 0) next

            q <- dbSendQuery(con,
                "SELECT departure_time FROM trip_updates WHERE date=? AND trip_id=? AND stop_id=?"
            )
            dbBind(q,
                list(
                    date,
                    trip,
                    stops[i]
                )
            )
            start_time <- dbFetch(q)$departure_time
            dbClearResult(q)
            if (length(start_time) == 0) next

            q <- dbSendQuery(con,
                "SELECT arrival_time FROM trip_updates WHERE date=? AND trip_id=? AND stop_id=?"
            )
            dbBind(q,
                list(
                    date,
                    trip,
                    stops[i+1]
                )
            )
            end_time <- dbFetch(q)$arrival_time
            dbClearResult(q)
            if (length(end_time) == 0) next

            if (end_time <= start_time) next

            q <- dbSendQuery(con,
                "INSERT INTO segments VALUES(?,?,?,?,?,?)"
            )
            dbBind(q,
                list(
                    date,
                    segment_ids[i],
                    trip,
                    start_time,
                    end_time,
                    end_time - start_time
                )
            )
            dbClearResult(q)
        }
    }
}

