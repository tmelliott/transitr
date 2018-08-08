## at some point, this will be a shiny app?

## but not yet 

##' View the realtime state of vehicles
##'
##' @title View Realtime Vehicle States
##' @param host the host name for the database
##' @return NULL
##' @author Tom Elliott
##' @export 
view_realtime <- function(host = Sys.getenv("transitr_host")) {
    vps <- get_vehicle_positions(host)
    
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



get_vehicle_positions <- function(host) {
    con <- RPostgreSQL::dbConnect(RPostgreSQL::PostgreSQL(), host = host, user="transitr", dbname="realtime")
    vps <- RPostgreSQL::dbReadTable(con, "vehicles")
    RPostgreSQL::dbDisconnect(con)

    vps
}
