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
    RSQLite::dbWriteTable(con, "road_segments", segments, append = TRUE)
    RSQLite::dbWriteTable(con, "intersections", intersections, append = TRUE)
    RSQLite::dbWriteTable(con, "shape_segments", sh_segs, append = TRUE)
    db_close(con)
}


merge_paths <- function(new, existing, n.thres = 3) {
    new.seg <- generate_segments(new, existing, n.thres)
    existing.seg <- generate_segments(existing, new, n.thres)

    NETWORK <- list(segments = tibble(id = integer(), from = integer(), 
                                      to = integer(), path = sf::st_sfc()),
                    intersections = tibble(id = integer(), position = sf::st_sfc(), 
                                           type = factor(character(), levels = c("int", "stop"))))

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


