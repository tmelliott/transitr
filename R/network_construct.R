construct_network <- function(nw, node_threshold = 0.01) {
    con <- db_connect(nw$database)

    nodes <- RSQLite::dbReadTable(con, "nodes")
    intersections <- RSQLite::dbReadTable(con, "intersections")
    segments <- RSQLite::dbReadTable(con, "road_segments")
    shape_nodes <- RSQLite::dbReadTable(con, "shape_nodes")

    ## Load STOPS -> nodes
    # checking each stop for an existing node within ~2m
    stops <- RSQLite::dbGetQuery(con,
        "SELECT * FROM stops"
    )
    if (interactive() && nrow(stops)) {
        cat("\n * Processing stops as nodes ...\n")
        pb <- utils::txtProgressBar(0, nrow(stops), style = 3)
    }
    for (i in seq_along(stops$stop_id)) {
        if (!is.na(stops$node_id[i])) next()
        if (interactive()) pb$up(pb$getVal() + 1)
        # NodeMatrix <- nodes[, c("node_lon", "node_lat")]
        # check for existing node
        node_id <- NA
        # if (nrow(NodeMatrix)) {
        #     # dists <- geosphere::distGeo(
        #     #     stops[i, c("stop_lon", "stop_lat")],
        #     #     NodeMatrix
        #     # )
        #     rad <- function(x) x*pi/180
        #     dists <- sqrt(
        #         ((rad(NodeMatrix[,1]) - rad(stops$stop_lon[i])) *
        #             cos(rad(stops$stop_lat[i])))^2 +
        #         (rad(NodeMatrix[,2]) - rad(stops$stop_lat[i]))^2
        #     ) * 6378137
        #     if (min(dists) < node_threshold)
        #         node_id <- nodes$node_id[which.min(dists)]
        # }

        stop_code <- stops$stop_code[i]
        matching_nodes <- stops$node_id[stops$stop_code == stop_code]
        if (length(matching_nodes))
            node_id <- matching_nodes[1]

        # if one doesn't exist, create it
        if (is.na(node_id)) {
            node_id <- max(nodes$node_id + 1L, 1L, na.rm = TRUE)
            nodes <- rbind(
                nodes,
                data.frame(
                    node_id = node_id,
                    node_type = 0L,
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
    res <- RSQLite::dbSendQuery(con, "DELETE FROM stops")
    RSQLite::dbClearResult(res)
    RSQLite::dbWriteTable(con, "stops", stops, append = TRUE)
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

    rownames(nodes) <- nodes$node_id
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
        rad <- function(x) x*pi/180
        shape$shape_dist_traveled <- cumsum(c(0, sqrt(
            ((rad(ShapeMat[-nrow(ShapeMat), 1]) - rad(ShapeMat[-1, 1])) *
                cos(rad(ShapeMat[-1, 2])))^2 +
            (rad(ShapeMat[-nrow(ShapeMat), 2]) - rad(ShapeMat[-1, 2]))^2
        ) * 6371000))

        ## New approach:
        # - use all nodes that are STOPS, and
        # - find all nodes that are INTERSECTIONS within ~5m of shape

        # There's an issue with CIRCULAR routes,
        # so need to actually fetch the stops and THEN the nodes ...

        qry <- RSQLite::dbSendQuery(con,
            paste(sep = "\n",
                "SELECT node_id FROM stop_times",
                "INNER JOIN stops ON stops.stop_id = stop_times.stop_id ",
                "WHERE trip_id=(",
                "    SELECT trip_id FROM trips WHERE shape_id=? LIMIT 1",
                ")",
                "ORDER BY stop_sequence"
            )
        )
        RSQLite::dbBind(qry, list(shape_id))
        stop_node_ids <- RSQLite::dbFetch(qry)$node_id
        RSQLite::dbClearResult(qry)

        stop_nodes <- nodes[as.character(stop_node_ids),]
        NodeMatrix <- stop_nodes[, c("node_lon", "node_lat")]
        # with(NodeMatrix, points(node_lon, node_lat))
        dmin <- 0
        node_dist <- apply(NodeMatrix, 1, function(n) {
            # d <- geosphere::distGeo(n, ShapeMat)
            # use equirectangular instead
            d <- sqrt(((rad(ShapeMat[,1]) - rad(n[1])) * cos(rad(n[2])))^2 +
                (rad(ShapeMat[,2]) - rad(n[2]))^2) * 6371000
            # - compute shape dist travelled of found nodes
            d[shape$shape_dist_traveled < dmin] <- NA
            dmin <<- shape$shape_dist_traveled[which.min(d)]
            dmin
        })
        if (any(is.na(node_dist))) {
            print(shape_id)
            print(stop_nodes)
            print(node_dist)
        }
        if (any(diff(node_dist) < 0)) {
            print("ERROR ERROR THERE'S AN ERROR")
            print(stop_nodes)
            print(node_dist)
        }
        ShapeNodes <- stop_nodes[!is.na(node_dist), ]
        node_dist <- node_dist[!is.na(node_dist)]

        # with(ShapeNodes, points(node_lon, node_lat, cex = 0.5, pch = 21, bg="red"))
        ShapeNodes$dist <- node_dist
        if (nrow(ShapeNodes) == 0) next
        #ShapeNodes <- ShapeNodes[order(ShapeNodes$dist),]
        ShapeNodes$seq <- seq_along(1:nrow(ShapeNodes))

        # with(ShapeNodes, points(node_lon, node_lat, cex = 0.5, pch = 21, bg="white"))
        # with(ShapeNodes, text(node_lon, node_lat, seq, col = "blue"))

        shape_nodes <- data.frame(
            shape_id = shape_id,
            node_id = ShapeNodes$node_id,
            node_sequence = ShapeNodes$seq,
            distance_traveled = ShapeNodes$dist
        )
        if (max(shape_nodes$distance_traveled) != max(shape$shape_dist_traveled)) {
            if (abs(max(shape_nodes$distance_traveled) - max(shape$shape_dist_traveled)) < 50)
                shape_nodes$distance_traveled[nrow(shape_nodes)] <-
                    max(shape$shape_dist_traveled)
            else {
                print("ERROR ERROR ERROR")
                print(utils::tail(shape))
                print(shape_nodes)
                next
            }
        }

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
        next_seg_id <- RSQLite::dbGetQuery(con,
            "SELECT MAX(road_segment_id)as ID FROM road_segments"
        )$ID + 1
        if (is.na(next_seg_id)) next_seg_id <- 1
        shape_segs$road_segment_id <-
            seq(next_seg_id, by = 1, length.out = nrow(shape_segs))
        RSQLite::dbWriteTable(con, "road_segments", shape_segs, append = TRUE)
    }
    if (exists("pb")) close(pb)

    db_close(con)
}
