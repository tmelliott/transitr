## at some point, this will be a shiny app?

## but not yet 

##' View the realtime state of vehicles
##'
##' @title View Realtime Vehicle States
##' @param host the host name for the database
##' @return NULL
##' @author Tom Elliott
##' @import shiny
##' @import ggplot2
##' @import leaflet
##' @export 
view_realtime <- function(host = Sys.getenv("transitr_host")) {
    vps <- get_vehicle_positions(host)

    options(transitr.host = host)
    app <- shinyApp(ui = rt_ui(), server = rt_server)
    runApp(app, launch.browser = FALSE)
}



get_vehicle_positions <- function(host) {
    con <- RPostgreSQL::dbConnect(RPostgreSQL::PostgreSQL(), 
                                  host = host, user = "transitr", 
                                  dbname = "realtime")
    vps <- RPostgreSQL::dbReadTable(con, "vehicles")
    RPostgreSQL::dbDisconnect(con)

    vps
}
