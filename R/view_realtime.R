## at some point, this will be a shiny app?

## but not yet 

##' View the realtime state of vehicles
##'
##' @title View Realtime Vehicle States
##' @param db SQLite database containing vehicle states, ETAs, etc
##' @return NULL
##' @author Tom Elliott
##' @export 
view_realtime <- function(db) {
    con <- db_connect(db)
    on.exit(db_close(con))

    vps <- RSQLite::dbReadTable(con, "vehicles")
    
    # p <- ggplot2::ggplot(NULL, ggplot2::aes(sort(vps$progress), 1:nrow(vps))) +
    #     ggplot2::geom_point() +
    #     ggplot2::xlim(0, 100)
    p <- ggplot2::ggplot(vps, ggplot2::aes(position_longitude, position_latitude, colour = speed*60*60/1000)) + 
        ggplot2::geom_point () +
        ggplot2::scale_color_viridis_c()
    dev.hold()
    print(p)
    dev.flush()

    invisible(NULL)
}
