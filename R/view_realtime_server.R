rt_server <- function(input, output, session) {
    # output$vehiclemap <- renderPlot(
    #     ggplot(get_vehicle_positions(options("transitr.host")), 
    #            aes(position_longitude, position_latitude, 
    #                colour = speed*60*60/1000)) + 
    #         geom_point () +
    #         scale_color_viridis_c()
    # )
    output$vehicles <- renderLeaflet({
        vps <- get_vehicle_positions(options("transitr.host"))
        m <- leaflet(vps)
        m <- addTiles(m)
        m <- addMarkers(m, ~position_longitude, ~position_latitude)
        m
    })
}

