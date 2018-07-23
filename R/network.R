##' Constructs a transit network from existing GTFS shapes.
##' Detects intersections and roads between them.
##'
##' @title Construct transit network
##' @param nw a \code{trgtfs} object
##' @return A \code{trgtfs} object, with additional network tables.
##' @author Tom Elliott
##' @export
construct <- function(nw) {
    create_network_tables(nw)
    tmp_construct_network(nw)
    
    invisible(nw)
}

load_shapes <- function(nw) {
    where <- ifelse(has_version_api(nw),
                    "WHERE `version`=(SELECT MAX(`version`) FROM `versions`)",
                    "")
    con <- db_connect(nw$database)
    shapes <- RSQLite::dbGetQuery(
        con,
        glue::glue(
            "SELECT `shape_id`, `shape_pt_lat`, `shape_pt_lon`, `shape_pt_sequence`
               FROM `shapes` {where}")
    )
    db_close(con)
    shapes$shape_id <- factor(shapes$shape_id)

    shapes_df_to_list(shapes)
}

create_network_tables <- function(nw) {
    con <- db_connect(nw$database)
    on.exit(db_close(con))
    if (RSQLite::dbExistsTable(con, "road_segments") ||
        RSQLite::dbExistsTable(con, "intersections") ||
        RSQLite::dbExistsTable(con, "trip_segments")) {
        stop("Database already has network tables")
    }
    
    res <- RSQLite::dbSendQuery(
        con,
        paste_nl(
            "CREATE TABLE road_segments (",
            "  road_segment_id INTEGER PRIMARY KEY,",
            "  int_from INTEGER,",
            "  int_to INTEGER,",
            "  length DOUBLE",
            ")"))
    RSQLite::dbClearResult(res)
        
    res <- RSQLite::dbSendQuery(
        con,
        paste_nl(
            "CREATE TABLE intersections (",
            "  intersection_id INTEGER PRIMARY KEY,",
            "  intersection_lat DOUBLE,",
            "  intersection_lon DOUBLE",
            ")"))
    RSQLite::dbClearResult(res)
        
    res <- RSQLite::dbSendQuery(
        con,
        paste_nl(
            "CREATE TABLE shape_segments (",
            "  shape_id TEXT,",
            "  road_segment_id INTEGER,",
            "  shape_road_sequence INTEGER,",
            "  distance_traveled DOUBLE",
            ")"))
    RSQLite::dbClearResult(res)
}

load_table <- function(nw, tbl) {
    con <- db_connect(nw)
    on.exit(db_close(con))
    RSQLite::dbReadTable(con, tbl)
}

load_road_segments <- function(nw) {
    load_table(nw, "road_segments")
}

load_intersections <- function(nw) {
    load_table(nw, "intersections")
}

load_shape_segments <- function(nw) {
    load_table(nw, "shape_segments")
}
