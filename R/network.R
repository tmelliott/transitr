##' Constructs a transit network from existing GTFS shapes.
##' Detects intersections and roads between them.
##'
##' @title Construct transit network
##' @param nw a \code{trgtfs} object
##' @return A \code{trgtfs} object, with additional network tables.
##' @author Tom Elliott
##' @export
construct <- function(nw) {
    shapes <- load_shapes(nw)
    
}

load_shapes <- function(nw) {
    where <- ifelse(has_version_api(nw),
                    "WHERE `version`=(SELECT MAX(`version`) FROM `versions`)",
                    "")
    shapes <- RSQLite::dbGetQuery(
        nw$connection,
        glue::glue(
            "SELECT `shape_id`, `shape_pt_sequence`, `shape_pt_lat`, `shape_pt_lon`
               FROM `shapes` {where}")
    )

    shapes_df_to_list(shapes)
}
