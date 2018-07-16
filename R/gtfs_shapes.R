create_shapes <- function(object) {
    con <- object$connection
    if (RSQLite::dbExistsTable(con, "shapes")) {
        stop("Shapes table already exists")
    }

    res <- dbSendQuery(
        con,
        paste_nl(
            "CREATE TABLE shapes (",
            "  shape_id TEXT,",
            "  shape_pt_lat DOUBLE,",
            "  shape_pt_lon DOUBLE,",
            "  shape_pt_sequence INTEGER,",
            "  shape_dist_traveled DOUBLE,",
            "  version DOUBLE,",
            "  PRIMARY KEY (shape_id, shape_pt_sequence)",
            ")"))
    dbClearResult(res)
}

update_shapes <- function(object, file) {
    existing <- dbGetQuery(object$connection,
                           "SELECT shape_id FROM shapes")
    tbl <- read.csv(file, header = TRUE)
    shapes <- data.frame(shape_id = as.character(tbl$shape_id),
                         shape_pt_lat = as.numeric(tbl$shape_pt_lat),
                         shape_pt_lon = as.numeric(tbl$shape_pt_lon),
                         shape_pt_sequence = as.integer(tbl$shape_pt_sequence),
                         shape_dist_traveled = as.double(tbl$shape_dist_traveled),
                         version = as.numeric(gsub(".+_v", "", tbl$shape_id)))
    shapes <- shapes[!shapes$shape_id %in% existing, ]
    dbWriteTable(object$connection, "shapes", shapes, append = TRUE)
}
