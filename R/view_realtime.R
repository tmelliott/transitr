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

    p <- ggplot2::ggplot(vps, ggplot2::aes(progress, vehicle_id)) +
        ggplot2::geom_point() +
        ggplot2::xlim(0, 100)
    dev.hold()
    print(p)
    dev.flush()

    invisible(NULL)
}
