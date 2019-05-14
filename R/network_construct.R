construct_network <- function(nw) {
    con <- db_connect(nw$database)

    intersections <- RSQLite::dbReadTable(con, "intersections")
    segments <- RSQLite::dbReadTable(con, "road_segments")
    shape_segments <- RSQLite::dbReadTable(con, "shape_segments")

    ## For each route ...
    routes <- RSQLite::dbGetQuery(con, "SELECT route_id FROM routes")

    if (interactive())
        pb <- utils::txtProgressBar(0, nrow(routes), style = 3)
    for (route in rev(routes$route_id)) {
        if (interactive()) pb$up(pb$getVal() + 1)
        
        trq <- RSQLite::dbSendQuery(con, "SELECT trip_id, shape_id FROM trips WHERE route_id=? LIMIT 1")
        RSQLite::dbBind(trq, list(route))
        res <- RSQLite::dbFetch(trq)
        RSQLite::dbClearResult(trq)

        if (nrow(shape_segments) > 0 && nrow(res) > 0 && res$shape_id %in% shape_segments$shape_id) next()

        ## ... get a trip and find stop sequence ...
        stq <- RSQLite::dbSendQuery(con, "SELECT stop_times.stop_id, stop_sequence, stop_lat, stop_lon FROM stop_times, stops WHERE trip_id=? AND stop_times.stop_id=stops.stop_id ORDER BY stop_sequence")
        RSQLite::dbBind(stq, list(res$trip_id))
        stops <- RSQLite::dbFetch(stq)
        RSQLite::dbClearResult(stq)

        ## ... and get the shape ...
        shq <- RSQLite::dbSendQuery(con, "SELECT shape_pt_lat, shape_pt_lon, shape_pt_sequence, shape_dist_traveled FROM shapes WHERE shape_id=? ORDER BY shape_pt_sequence")
        RSQLite::dbBind(shq, list(res$shape_id))
        shape <- RSQLite::dbFetch(shq)
        RSQLite::dbClearResult(shq)

        ## ... then create "intersections" at stops ... 
        intid = NA
        if (nrow(stops) == 0) next()
        segs <- integer(nrow(stops) - 1L)
        for (i in 1:nrow(stops)) {
            previd <- intid
            si <- stops[i, ]
            inti <- intersections[intersections$intersection_lat == si$stop_lat & 
                                  intersections$intersection_lon == si$stop_lon, ]
            if (nrow(inti) == 0) {
                ## create new intesection
                intid <- max(c(0, intersections$intersection_id)) + 1
                intersections <- rbind(intersections, 
                    data.frame(intersection_id = intid, intersection_lat = si$stop_lat, intersection_lon = si$stop_lon))
            } else {
                ## use existing
                intid <- inti$intersection_id[1]
            }

            if (i > 1) {
                ## ... and segments between intersections.
                segi <- segments[segments$int_from == previd & segments$int_to == intid, ]
                if (nrow(segi) == 0) {
                    ## create new segment
                    segid <- max(c(0, segments$road_segment_id)) + 1
                    segments <- rbind(segments,
                        data.frame(road_segment_id = segid, int_from = previd, int_to = intid, length = 0))
                } else {
                    ## use existing
                    segid <- segi$road_segment_id
                }
                segs[i-1] <- segid
            }
        }

        ## Finally, create shape_segments - which requires distances
        if (any(is.na(shape$shape_dist_traveled))) {
            if (nrow(shape) == 1) {
                shape$shape_dist_traveled <- 0
            } else if (nrow(shape) == 2) {
                shape$shape_dist_traveled <- c(
                    0,
                    geosphere::distGeo(
                        shape[1, 2L:1L, drop = FALSE], 
                        shape[2, 2L:1L, drop = FALSE]
                    )
                )
            } else {
                shape$shape_dist_traveled <- c(0, cumsum(geosphere::distGeo(shape[,2L:1L])))
            }
        }

        siprev <- 0
        si <- 1
        for (i in 1:length(segs)) {
            ## first we need to find the closest point in shape to stop ...
            segi <- segments[segments$road_segment_id == segs[i], ]
            iprev <- intersections[intersections$intersection_id == segi$int_from, ]
            icur <- intersections[intersections$intersection_id == segi$int_to, ]

            siprev <- si + which.min(geosphere::distGeo(iprev[, 3L:2L], shape[si:nrow(shape), 2L:1L])) - 1
            si <- siprev + which.min(geosphere::distGeo(icur[, 3L:2L], shape[siprev:nrow(shape), 2L:1L])) - 1
            ## ... then use shape_dist_traveled to insert segment
            shape_segments <- rbind(shape_segments, 
                data.frame(shape_id = res$shape_id,
                           road_segment_id = segs[i],
                           shape_road_sequence = i,
                           distance_traveled = shape$shape_dist_traveled[siprev]))
            if (segi$length == 0)
                segments[segments$road_segment_id == segs[i], "length"] <- 
                    shape$shape_dist_traveled[si] - shape$shape_dist_traveled[siprev]
        }
    }
    if (interactive()) close(pb)

    RSQLite::dbWriteTable(con, "road_segments", segments, overwrite = TRUE)
    RSQLite::dbWriteTable(con, "intersections", intersections, overwrite = TRUE)
    RSQLite::dbWriteTable(con, "shape_segments", shape_segments, overwrite = TRUE)

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
