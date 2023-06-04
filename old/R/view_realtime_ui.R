rt_ui <- function() {
    fluidPage(
        tags$head(
            tags$style(HTML("
                html, body {
                    height: 100%
                }
            "))
        ),
        style = "height: 100%",
        fluidRow(
            style = "height: 100%",
            column(2,
                textOutput("count")
            ),
            column(10,
                style = "height: 100%",
                leafletOutput("vehicles", height = "100%")
            )
        )
    )
}
