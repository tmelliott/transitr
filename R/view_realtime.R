## at some point, this will be a shiny app?

## but not yet 

#' @import ggplot2
#' @export
view_realtime <- function(db) {
    con <- db_connect(db)
    on.exit(db_close(con))

    vps <- RSQLite::dbReadTable(con, "vehicles")

    p <- ggplot(vps, aes(progress, vehicle_id)) +
        geom_point() +
        xlim(0, 100)
    dev.hold()
    print(p)
    dev.flush()

    invisible(NULL)
}
