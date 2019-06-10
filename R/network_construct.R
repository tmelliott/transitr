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


merge_paths <- function(new, existing, n.thres = 3) {
    # new.seg <- generate_segments(new, existing, n.thres)
    # existing.seg <- generate_segments(existing, new, n.thres)

    # NETWORK <- list(segments = tibble(id = integer(), from = integer(), 
    #                                   to = integer(), path = sf::st_sfc()),
    #                 intersections = tibble(id = integer(), position = sf::st_sfc(), 
    #                                        type = factor(character(), levels = c("int", "stop"))))

    # new.seg.geom <- do.call(st_sfc, tapply(seq_along(1:nrow(new.seg)), new.seg[, 5],
    #                                        function(i) new.seg[i, 2:3] %>% st_linestring))
    # existing.seg.geom <- do.call(st_sfc, tapply(seq_along(1:nrow(existing.seg)), existing.seg[, 5],
    #                                             function(i) existing.seg[i, 2:3] %>% st_linestring))

    # plot(new, type = "l")
    # lines(new.seg[, 2:3])
    # plot(new.seg.geom[[2]])
    # plot(existing.seg.geom)

    # NZ: 5759
    CEN <- c(36.8485, 174.7633)

    nl <- st_linestring(new %>% as.matrix)
    el <- st_linestring(existing %>% as.matrix)
    bl <- st_sfc(nl, el, crs = 4326)
    blt <- st_transform(bl, sprintf("+proj=eqc +lat_0=%s +lon_0=%s", CEN[1], CEN[2]))
    plot(blt)

    blb <- st_buffer(blt, 20)

    plot(blb)
    plot(st_difference(blb[[1]], blb[[2]]), col=  'red', add=TRUE)
    plot(st_difference(blb[[2]], blb[[1]]), col= 'blue', add=TRUE)
    plot(st_intersection(blb[[1]], blb[[2]]), col = "yellow", add = TRUE)

    int <- st_intersection(blb[[1]], blb[[2]])

    ## take EXISTING and split it up
    p <- ggplot(existing, aes(shape_pt_lon, shape_pt_lat)) + geom_path(color = "gray") +
        geom_path(data = new, color = "gray") + coord_fixed(ratio = 1.6)
    for (si in unique(existing.seg[existing.seg[,4] == 1, 5])) {
        wi <- which(existing.seg[,5] == si)
        ei <- existing[wi,]
        for (xi in unique(new.seg[new.seg[,4] == 1, 5])) {
            vi <- which(new.seg[,5] == xi)
            ni <- new[vi,]
            xr <- extendrange(ei$shape_pt_lon)
            yr <- extendrange(ei$shape_pt_lat)
            px <- p + xlim(xr[1], xr[2]) + ylim(yr[1], yr[2]) +
                geom_path(data = ei, colour = "orangered", lwd = 2) +
                geom_path(data = ni, color = "blue")
            geosphere::distGeo(apply(ni, 2, median), apply(ei, 2, median))
            if (geosphere::distGeo(ni[1,], ei[1,]) < 50 &
                geosphere::distGeo(ni[nrow(ni),], ei[nrow(ei),]) < 50) {
                break
            }
        }

        START <- END <- NULL

        ## START INTERSECTION
        if (si > 1) {
            if (geosphere::distGeo(ei[1,], ni[1,]) < 1) {
                START <- ei[1,]
            } else {
                px <- p + xlim(xr[1], xr[2]) + ylim(yr[1], yr[2]) +
                    geom_path(data = ei, colour = "orangered", lwd = 2) +
                    geom_path(data = ni, color = "blue") +
                    geom_point(data = existing[min(wi)+c(-1,0),], color = "black") +
                    geom_point(data = new[min(vi)+c(-1,0),], color = "blue")
                px
                ## which is closer to the end of the previous segments?
                dee <- geosphere::distGeo(ei[1,], existing[min(wi)-1,])
                den <- geosphere::distGeo(ei[1,], new[min(vi)-1,])
                dne <- geosphere::distGeo(ni[1,], existing[min(wi)-1,])
                dnn <- geosphere::distGeo(ni[1,], new[min(vi)-1,])
                if (dee + den < dne + dnn) {
                    START <- existing.seg[min(wi),2:3]
                } else {
                    START <- new.seg[min(vi),2:3]
                }
                START <- rbind(START) %>% as.tibble
                colnames(START) <- c("shape_pt_lon", "shape_pt_lat")
            }

            xr <- extendrange(c(ei$shape_pt_lon[1:5], ni$shape_pt_lon[1:5]), f=0.5)
            yr <- extendrange(c(ei$shape_pt_lat[1:5], ni$shape_pt_lat[1:5]), f=0.5)
            px <- p + xlim(xr[1], xr[2]) + ylim(yr[1], yr[2]) +
                geom_path(data = ei, colour = "orangered", lwd = 2) +
                geom_path(data = ni, color = "blue") +
                geom_point(data = START)
            print(px)
            grid::grid.locator()
        } else {
            ## it's the start of a route (i.e., a stop) AND overlaps a segment ... 
            #... it's a PARTIAL SEGMENT
        }

        ## END INTERSECTION
        if (si < max(existing.seg[,5])) {
            if (geosphere::distGeo(ei[nrow(ei),], ni[nrow(ni),]) < 1) {
                END <- ei[nrow(ei),]
            } else {
                px <- p + xlim(xr[1], xr[2]) + ylim(yr[1], yr[2]) +
                    geom_path(data = ei, colour = "orangered", lwd = 2) +
                    geom_path(data = ni, color = "blue") +
                    geom_point(data = existing[max(wi)+c(1,0),], color = "black") +
                    geom_point(data = new[max(vi)+c(1,0),], color = "blue")
                px
                ## which is closer to the end of the previous segments?
                dee <- geosphere::distGeo(ei[nrow(ei),], existing[max(wi)+1,])
                den <- geosphere::distGeo(ei[nrow(ei),], new[max(vi)+1,])
                dne <- geosphere::distGeo(ni[nrow(ni),], existing[max(wi)+1,])
                dnn <- geosphere::distGeo(ni[nrow(ni),], new[max(vi)+1,])
                if (dee + den < dne + dnn) {
                    END <- existing.seg[max(wi),2:3]
                } else {
                    END <- new.seg[max(vi),2:3]
                }
                END <- rbind(END) %>% as.tibble
                colnames(END) <- c("shape_pt_lon", "shape_pt_lat")
            }
            xr <- extendrange(c(tail(ei$shape_pt_lon, 5), tail(ni$shape_pt_lon, 5)), f=0.5)
            yr <- extendrange(c(tail(ei$shape_pt_lat, 5), tail(ni$shape_pt_lat, 5)), f=0.5)
            px <- p + xlim(xr[1], xr[2]) + ylim(yr[1], yr[2]) +
                geom_path(data = ei, colour = "orangered", lwd = 2) +
                geom_path(data = ni, color = "blue") +
                geom_point(data = END)
            print(px)
            grid::grid.locator()
        } else {
            ## it's at the end of a route (i.e., a stop) AND overlaps a segment ...
            #... it's a PARTIAL SEGMENT
        }

        xr <- extendrange(ei$shape_pt_lon)
        yr <- extendrange(ei$shape_pt_lat)
        px <- p + xlim(xr[1], xr[2]) + ylim(yr[1], yr[2]) +
            geom_path(data = ei, colour = "orangered", lwd = 2) +
            geom_path(data = ni, color = "blue") +
            geom_point(data = rbind(START, END))
        print(px)

        grid::grid.locator()
    }

    ## first, overlapping segments
    

    ## then, the remainder


    ## then create remaining segments for NEW
    new.segments <- tapply(seq_along(segi), segi, function(i) new[i, ])

    ## then for each segment, break EXISTING up
    

    segi
}

generate_segments <- function(target, ref, n.thres, unique.only = FALSE) {
    ## calculate distance from point in new to nearest point in existing
    dsx <- apply(target, 1, geosphere::dist2Line, line = ref)

    ## remove any short segments
    ds <- dsx[1,]
    di <- ds < 3
    d1 <- d2 <- integer(length(ds))
    for (i in 2:length(d1)) if (di[i] == di[i-1]) d1[i] <- d1[i-1] + 1
    for (i in (length(d2)-1):1) if (di[i] == di[i+1]) d2[i] <- d2[i+1] + 1

    if (any(d1 + d2 < n.thres)) {
        di[d1 + d2 < n.thres] <- !di[d1 + d2 < n.thres]
        d1 <- d2 <- integer(length(ds))
        for (i in 2:length(d1)) if (di[i] == di[i-1]) d1[i] <- d1[i-1] + 1
        for (i in (length(d2)-1):1) if (di[i] == di[i+1]) d2[i] <- d2[i+1] + 1    
    }

    segii <- cumsum(d1 == 0)
    return(cbind(t(dsx), di, segii))


    if (unique.only) 
        d1[di] <- NA

    ## figure out segment indexes
    segi <- d1
    segi[!is.na(d1)] <- cumsum(d1[!is.na(d1)] == 0)

    ## add end points
    segments <- vector("list", length(unique(segi)))
    intersection <- NULL
    for (seg in unique(segi)) {
        w <- which(segi == seg)
        this.seg <- rbind(intersection, target[w, ])
        if (seg < max(segi)) {
            w1 <- max(w)
            s0 <- target[w0,]
            
            if (ds[w1] < ds[w1+1]) {
                ## find point at intersection and insert it
                intersection <- dsx[2:3, w1]
            } else {
                intersection <- dsx[2:3, w1+1]
            }
            this.seg <- rbind(this.seg, intersection)
        }
        segments[[seg]] <- cbind(this.seg, i = seg)
    }

    segs <- lapply(segments, as.tibble) %>% bind_rows
    attr(segs, "index") <- tapply(segi, segii, function(x) if (any(is.na(x))) NA else x[1])
    segs
}

delete_stops <- function(path, stops) {
    d <- apply(path[, c("shape_pt_lon", "shape_pt_lat")] %>% as.matrix, 1, 
        function(x) min(geosphere::distGeo(x, p2 = stops)))
    for (i in 1:(length(d)-1)) if (d[i+1] < 2) d[i] <- NA
    path <- path[!is.na(d), ]
    d <- d[!is.na(d)]
    # interpolate 
    for (i in which(d == 0)) {
        if (i == 1 || i == length(d)) next
        path[i, c("shape_pt_lon", "shape_pt_lat")] <- 
            geosphere::dist2Line(path[i, c("shape_pt_lon", "shape_pt_lat")], 
                path[c(i-1, i+1), c("shape_pt_lon", "shape_pt_lat")] %>% as.matrix)[, 2:3]
    }
    path
}


