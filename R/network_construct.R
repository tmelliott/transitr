construct_network <- function(nw) {
    con <- db_connect(nw$database)

    intersections <- dbReadTable(con, "intersections")
    segments <- dbReadTable(con, "road_segments")
    shape_segments <- dbReadTable(con, "shape_segments")

    ## For each route ...
    routes <- dbGetQuery(con, "SELECT route_id FROM routes")
    trq <- dbSendQuery(con, "SELECT trip_id, shape_id FROM trips WHERE route_id=? LIMIT 1")
    stq <- dbSendQuery(con, "SELECT stop_times.stop_id, stop_sequence, stop_lat, stop_lon FROM stop_times, stops WHERE trip_id=? AND stop_times.stop_id=stops.stop_id ORDER BY stop_sequence")
    shq <- dbSendQuery(con, "SELECT shape_pt_lat, shape_pt_lon, shape_pt_sequence FROM shapes WHERE shape_id=? ORDER BY shape_pt_sequence")

    pb <- txtProgressBar(0, nrow(routes), style = 3)
    for (route in routes$route_id) {
        pb$up(pb$getVal() + 1)
        dbBind(trq, list(route))
        res <- dbFetch(trq)
        dbClearResult(trq)

        ## ... get a trip and find stop sequence ...
        dbBind(stq, list(res$trip_id))
        stops <- dbFetch(stq)
        dbClearResult(stq)

        ## ... and get the shape ...
        dbBind(shq, list(res$shape_id))
        shape <- dbFetch(shq)
        dbClearResult(shq)

        ## ... then create "intersections" at stops ... 
        intid = NA
        segs <- integer(nrow(stops) - 1)
        for (i in 1:nrow(stops)) {
            previd <- intid
            si <- stops[i, ]
            inti <- intersections[intersections$intersection_lat == si$stop_lat & 
                                  intersections$intersection_lon == si$stop_lon, ]
            if (nrow(inti) == 0) {
                ## create new intesection
                intid <- max(1, max(intersections$intersection_id) + 1)
                intersections <- rbind(intersections, 
                    data.frame(intersection_id = max(intersections$intersection_id)+1,
                               intersection_lat = si$stop_lat,
                               intersection_lon = si$stop_lon))
            } else {
                ## use existing
                intid <- inti$intersection_id[1]
            }

            if (i > 1) {
                ## ... and segments between intersections.
                segi <- segments[segments$int_from == previd & segments$int_to == intid, ]
                if (nrow(segi) == 0) {
                    ## create new segment
                    segid <- max(1, max(segments$road_segment_id) + 1)
                    segments <- rbind(segments,
                        data.frame(road_segment_id = segid, int_from = previd, int_to = intid, length = 0))
                } else {
                    ## use existing
                    segid <- segi$road_segment_id
                }
                segs[i-1] <- segid
            }
        }

        ## Finally, create shape_segments

    }

    db_close(con)
}

tmp_construct_network <- function(nw) {
    shapes <- load_shapes(nw)

    con <- db_connect(nw$database)
    segments <- RSQLite::dbReadTable(con, "road_segments")
    intersections <- RSQLite::dbReadTable(con, "intersections")
    sh_segs <- RSQLite::dbReadTable(con, "shape_segments")
    db_close(con)

    ## To-do: check for duplicate intersections
    for (shape in shapes) {
        ni <- nrow(intersections)
        x <- shape[c(1, nrow(shape)), ]
        x <- data.frame(intersection_id = ni + 1:2, intersection_lon = x[,1], intersection_lat = x[,2])
        intersections <- rbind(intersections, x)
        segments <- rbind(segments, 
            data.frame(road_segment_id = nrow(segments) + 1, int_from = ni + 1, int_to = ni + 2, 
                       length = numeric(1)))
        sh_segs <- rbind(sh_segs, 
            data.frame(shape_id = attr(shape, "id"), road_segment_id = nrow(segments) + 1, 
                       shape_road_sequence = 1, distance_traveled= 0))
    }

    con <- db_connect(nw$database)
    RSQLite::dbWriteTable(con, "road_segments", segments, overwrite = TRUE)
    RSQLite::dbWriteTable(con, "intersections", intersections, overwrite = TRUE)
    RSQLite::dbWriteTable(con, "shape_segments", sh_segs, overwrite = TRUE)
    db_close(con)
}
