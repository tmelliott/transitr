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


merge_paths <- function(new, existing) {
    ## calculate distance from point in new to nearest point in existing
    ds <- apply(new, 1, geosphere::dist2Line, line = existing)
    ds[1, ]
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
