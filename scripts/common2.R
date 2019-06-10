get_segment_data <- function(routes) {
    con <- dbConnect(SQLite(), "fulldata.db")
    on.exit(dbDisconnect(con))
    if (missing(routes)) {
        segments <- con %>% tbl("road_segments")
    } else {
        rids <- con %>% tbl("routes") %>%
            filter(route_short_name %in% routes &
                (route_long_name %like% "%To City%" |
                 route_long_name %like% "%To Britomart%" |
                 route_long_name %like% "%To Mayoral%" |
                 route_long_name %like% "%To Auckland Universities%" |
                 route_long_name %like% "%To Hibiscus Coast Station")) %>%
            collect %>% pluck("route_id") %>% unique
        sids <- con %>% tbl("trips") %>%
            filter(route_id %in% rids) %>%
            collect %>% pluck("shape_id") %>% unique
        segids <- con %>% tbl("shape_segments") %>% 
            filter(shape_id %in% sids) %>%
            collect %>% pluck("road_segment_id") %>% unique
        segments <- con %>% tbl("road_segments") %>%
            filter(road_segment_id %in% segids)
    }
    intersections <- con %>% tbl("intersections")
    segments <- segments %>% 
        inner_join(intersections, by = c("int_from" = "intersection_id"), suffix = c("", "_start")) %>%
        inner_join(intersections, by = c("int_to" = "intersection_id"), suffix = c("", "_end")) %>%
        select(road_segment_id, length, intersection_lat, intersection_lon, intersection_lat_end, intersection_lon_end) %>%
        collect
    segments
}
