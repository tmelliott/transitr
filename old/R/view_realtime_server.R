rt_server <- function(input, output, session) {
    # output$vehiclemap <- renderPlot(
    #     ggplot(get_vehicle_positions(options("transitr.host")), 
    #            aes(position_longitude, position_latitude, 
    #                colour = speed*60*60/1000)) + 
    #         geom_point () +
    #         scale_color_viridis_c()
    # )
    vps <- get_vehicle_positions(options("transitr.host"))
    vps <- vps[Sys.time() - vps$timestamp < 120, ]
    output$count <- renderText(sprintf("Showing %i vehicle locations", nrow(vps)))
    output$vehicles <- renderLeaflet({
        m <- leaflet(vps)
        m <- addTiles(m)
        m <- addMarkers(m, ~position_longitude, ~position_latitude)
        m
    })
}

