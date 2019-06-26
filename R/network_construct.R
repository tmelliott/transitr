construct_network <- function(nw, node_threshold = 0.5) {
    con <- db_connect(nw$database)

    nodes <- RSQLite::dbReadTable(con, "nodes")
    intersections <- RSQLite::dbReadTable(con, "intersections")
    segments <- RSQLite::dbReadTable(con, "road_segments")
    shape_nodes <- RSQLite::dbReadTable(con, "shape_nodes")

    ## Load STOPS -> nodes
    # checking each stop for an existing node within ~2m
    stops <- RSQLite::dbGetQuery(con,
        "SELECT * FROM stops WHERE node_id is null"
    )
    if (interactive() && nrow(stops)) {
        cat("\n * Processing stops as nodes ...\n")
        pb <- utils::txtProgressBar(0, nrow(stops), style = 3)
    }
    for (i in seq_along(stops$stop_id)) {
        if (interactive()) pb$up(pb$getVal() + 1)
        NodeMatrix <- nodes[, c("node_lon", "node_lat")]
        # check for existing node
        node_id <- NA
        if (nrow(NodeMatrix)) {
            # dists <- geosphere::distGeo(
            #     stops[i, c("stop_lon", "stop_lat")],
            #     NodeMatrix
            # )
            rad <- function(x) x*pi/180
            dists <- sqrt(
                ((rad(NodeMatrix[,1]) - rad(stops$stop_lon[i])) * 
                    cos(rad(stops$stop_lat[i])))^2 +
                (rad(NodeMatrix[,2]) - rad(stops$stop_lat[i]))^2
            ) * 6378137
            if (min(dists) < node_threshold)
                node_id <- nodes$node_id[which.min(dists)]
        }
        # if one doesn't exist, create it
        if (is.na(node_id)) {
            node_id <- max(nodes$node_id + 1L, 1L)
            nodes <- rbind(
                nodes,
                data.frame(
                    node_id = node_id,
                    node_type = 0,
                    node_lon = stops$stop_lon[i],
                    node_lat = stops$stop_lat[i]
                )
            )
        }
        
        # then add node ID to stop
        stops$node_id[i] <- node_id
    }
    if (exists("pb")) close(pb)
    RSQLite::dbWriteTable(con, "nodes", nodes, overwrite = TRUE)
    RSQLite::dbWriteTable(con, "stops", stops, overwrite = TRUE)
    rm("stops")

    ## For each SHAPE that doesn't exist in shape_segments
    shapes <- RSQLite::dbGetQuery(con,
        "SELECT DISTINCT shape_id FROM shapes WHERE shape_id NOT IN (SELECT DISTINCT shape_id FROM shape_nodes)"
    )$shape_id
    NodeMatrix <- nodes[, c("node_lon", "node_lat")]
    if (interactive() && length(shapes)) {
        cat("\n * Processing shapes ...\n")
        pb <- utils::txtProgressBar(0, length(shapes), style = 3)
    }
    for (shape_id in shapes) {
        if (interactive()) pb$up(pb$getVal() + 1)
        # - load
        qry <- RSQLite::dbSendQuery(con, "SELECT * FROM shapes WHERE shape_id=?")
        RSQLite::dbBind(qry, list(shape_id))
        shape <- RSQLite::dbFetch(qry)
        RSQLite::dbClearResult(qry)

        ShapeMat <- as.matrix(shape[, c("shape_pt_lon", "shape_pt_lat")])
        # plot(ShapeMat, type="l")
        # points(NodeMatrix, pch = 21, cex = 0.5, col = "black", bg="white")

        # shape distance
        # distGeo now adds an NA
        shape$shape_dist_traveled <- 
            c(0, cumsum(geosphere::distGeo(ShapeMat)[-nrow(ShapeMat)]))

        
        # - find all NODES that are within ~5m of shape
        node_dist <- apply(NodeMatrix, 1, function(n) {
            # d <- geosphere::distGeo(n, ShapeMat)
            # use equirectangular instead
            rad <- function(x) x*pi/180
            d <- sqrt(((rad(ShapeMat[,1]) - rad(n[1])) * cos(rad(n[2])))^2 +
                (rad(ShapeMat[,2]) - rad(n[2]))^2) * 6378137
            if (min(d) > 5) return(NA)
            # - compute shape dist travelled of found nodes
            shape$shape_dist_traveled[which.min(d)]
        })
        ShapeNodes <- nodes[!is.na(node_dist), ]
        node_dist <- node_dist[!is.na(node_dist)]

        # with(ShapeNodes, points(node_lon, node_lat, cex = 0.5, pch = 21, bg="red"))

        ## Filter out unscheduled stops
        qry <- RSQLite::dbSendQuery(con, 
            paste_nl(
                "SELECT node_id FROM stops WHERE stop_id IN ",
                "  (SELECT stop_id FROM stop_times",
                "  WHERE trip_id=(SELECT trip_id FROM trips WHERE shape_id=? LIMIT 1))"
            )
        )
        RSQLite::dbBind(qry, list(shape_id))
        stopnodes <- RSQLite::dbFetch(qry)$node_id
        RSQLite::dbClearResult(qry)
        ShapeNodes$dist <- node_dist
        ShapeNodes <- ShapeNodes[ShapeNodes$node_id %in% stopnodes | ShapeNodes$node_type == 1,]
        if (nrow(ShapeNodes) == 0) next
        ShapeNodes <- ShapeNodes[order(ShapeNodes$dist),]
        ShapeNodes$seq <- seq_along(1:nrow(ShapeNodes))

        # with(ShapeNodes, points(node_lon, node_lat, cex = 0.5, pch = 21, bg="white"))
        # with(ShapeNodes, text(node_lon, node_lat, seq, col = "blue"))
        
        shape_nodes <- data.frame(
            shape_id = shape_id,
            node_id = ShapeNodes$node_id,
            node_sequence = ShapeNodes$seq,
            distance_traveled = ShapeNodes$dist
        )

        # - insert into shape_nodes
        RSQLite::dbWriteTable(con, "shape_nodes", shape_nodes, append = TRUE)
        
        # - find segments between found nodes
        shape_segs <- data.frame(
            road_segment_id = NA,
            node_from = shape_nodes$node_id[-nrow(shape_nodes)],
            node_to = shape_nodes$node_id[-1],
            length = diff(shape_nodes$distance_traveled)
        )
        qry <- RSQLite::dbSendQuery(con,
            "SELECT road_segment_id FROM road_segments WHERE node_from=? AND node_to=?"
        )
        for (i in seq_along(1:nrow(shape_segs))) {
            RSQLite::dbBind(qry, 
                list(shape_segs$node_from[i], shape_segs$node_to[i])
            )
            seg <- RSQLite::dbFetch(qry)$road_segment_id
            if (length(seg)) shape_segs$road_segment_id[i] <- seg            
        }
        RSQLite::dbClearResult(qry)

        # insert new segments
        shape_segs <- shape_segs[is.na(shape_segs$road_segment_id), ]
        RSQLite::dbWriteTable(con, "road_segments", shape_segs, append = TRUE)
    }
    if (exists("pb")) close(pb)


    # routes <- RSQLite::dbGetQuery(con, "SELECT route_id FROM routes")
    # if (interactive())
    #     pb <- utils::txtProgressBar(0, nrow(routes), style = 3)
    # for (route in rev(routes$route_id)) {
    #     if (interactive()) pb$up(pb$getVal() + 1)
        
    #     trq <- RSQLite::dbSendQuery(con, "SELECT trip_id, shape_id FROM trips WHERE route_id=? LIMIT 1")
    #     RSQLite::dbBind(trq, list(route))
    #     res <- RSQLite::dbFetch(trq)
    #     RSQLite::dbClearResult(trq)

    #     if (nrow(shape_segments) > 0 && nrow(res) > 0 && res$shape_id %in% shape_segments$shape_id) next()

    #     ## ... get a trip and find stop sequence ...
    #     stq <- RSQLite::dbSendQuery(con, "SELECT stop_times.stop_id, stop_sequence, stop_lat, stop_lon FROM stop_times, stops WHERE trip_id=? AND stop_times.stop_id=stops.stop_id ORDER BY stop_sequence")
    #     RSQLite::dbBind(stq, list(res$trip_id))
    #     stops <- RSQLite::dbFetch(stq)
    #     RSQLite::dbClearResult(stq)

    #     ## ... and get the shape ...
    #     shq <- RSQLite::dbSendQuery(con, "SELECT shape_pt_lat, shape_pt_lon, shape_pt_sequence, shape_dist_traveled FROM shapes WHERE shape_id=? ORDER BY shape_pt_sequence")
    #     RSQLite::dbBind(shq, list(res$shape_id))
    #     shape <- RSQLite::dbFetch(shq)
    #     RSQLite::dbClearResult(shq)

    #     ## ... then create "intersections" at stops ... 
    #     intid = NA
    #     if (nrow(stops) == 0) next()
    #     segs <- integer(nrow(stops) - 1L)
    #     for (i in 1:nrow(stops)) {
    #         previd <- intid
    #         si <- stops[i, ]
    #         inti <- intersections[intersections$intersection_lat == si$stop_lat & 
    #                               intersections$intersection_lon == si$stop_lon, ]
    #         if (nrow(inti) == 0) {
    #             ## create new intesection
    #             intid <- max(c(0, intersections$intersection_id)) + 1
    #             intersections <- rbind(intersections, 
    #                 data.frame(intersection_id = intid, intersection_lat = si$stop_lat, intersection_lon = si$stop_lon))
    #         } else {
    #             ## use existing
    #             intid <- inti$intersection_id[1]
    #         }

    #         if (i > 1) {
    #             ## ... and segments between intersections.
    #             segi <- segments[segments$int_from == previd & segments$int_to == intid, ]
    #             if (nrow(segi) == 0) {
    #                 ## create new segment
    #                 segid <- max(c(0, segments$road_segment_id)) + 1
    #                 segments <- rbind(segments,
    #                     data.frame(road_segment_id = segid, int_from = previd, int_to = intid, length = 0))
    #             } else {
    #                 ## use existing
    #                 segid <- segi$road_segment_id
    #             }
    #             segs[i-1] <- segid
    #         }
    #     }

    #     ## Finally, create shape_segments - which requires distances
    #     if (any(is.na(shape$shape_dist_traveled))) {
    #         if (nrow(shape) == 1) {
    #             shape$shape_dist_traveled <- 0
    #         } else if (nrow(shape) == 2) {
    #             shape$shape_dist_traveled <- c(
    #                 0,
    #                 geosphere::distGeo(
    #                     shape[1, 2L:1L, drop = FALSE], 
    #                     shape[2, 2L:1L, drop = FALSE]
    #                 )
    #             )
    #         } else {
    #             shape$shape_dist_traveled <- c(0, cumsum(geosphere::distGeo(shape[,2L:1L])))
    #         }
    #     }

    #     siprev <- 0
    #     si <- 1
    #     for (i in 1:length(segs)) {
    #         ## first we need to find the closest point in shape to stop ...
    #         segi <- segments[segments$road_segment_id == segs[i], ]
    #         iprev <- intersections[intersections$intersection_id == segi$int_from, ]
    #         icur <- intersections[intersections$intersection_id == segi$int_to, ]

    #         siprev <- si + which.min(geosphere::distGeo(iprev[, 3L:2L], shape[si:nrow(shape), 2L:1L])) - 1
    #         si <- siprev + which.min(geosphere::distGeo(icur[, 3L:2L], shape[siprev:nrow(shape), 2L:1L])) - 1
    #         ## ... then use shape_dist_traveled to insert segment
    #         shape_segments <- rbind(shape_segments, 
    #             data.frame(shape_id = res$shape_id,
    #                        road_segment_id = segs[i],
    #                        shape_road_sequence = i,
    #                        distance_traveled = shape$shape_dist_traveled[siprev]))
    #         if (segi$length == 0)
    #             segments[segments$road_segment_id == segs[i], "length"] <- 
    #                 shape$shape_dist_traveled[si] - shape$shape_dist_traveled[siprev]
    #     }
    # }
    # if (interactive()) close(pb)

    # RSQLite::dbWriteTable(con, "road_segments", segments, overwrite = TRUE)
    # RSQLite::dbWriteTable(con, "intersections", intersections, overwrite = TRUE)
    # RSQLite::dbWriteTable(con, "shape_segments", shape_segments, overwrite = TRUE)

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
