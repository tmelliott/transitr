create_shapes <- function(con) {
    if (RSQLite::dbExistsTable(con, "shapes")) {
        stop("Shapes table already exists")
    }

    res <- RSQLite::dbSendQuery(
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
    RSQLite::dbClearResult(res)
}

update_shapes <- function(object, file) {
    existing <- RSQLite::dbGetQuery(object$connection,
                           "SELECT shape_id FROM shapes")
    tbl <- read.csv(file, header = TRUE)
    shapes <- data.frame(shape_id = as.character(tbl$shape_id),
                         shape_pt_lat = as.numeric(tbl$shape_pt_lat),
                         shape_pt_lon = as.numeric(tbl$shape_pt_lon),
                         shape_pt_sequence = as.integer(tbl$shape_pt_sequence),
                         shape_dist_traveled = as.double(tbl$shape_dist_traveled),
                         version = as.numeric(gsub(".+_v", "", tbl$shape_id)))
    shapes <- shapes[!shapes$shape_id %in% existing, ]
    RSQLite::dbWriteTable(object$connection, "shapes", shapes, append = TRUE)
}

check_shapes <- function(con) {
    identical(RSQLite::dbGetQuery(con, "PRAGMA table_info(shapes)"),
              data.frame(cid = 0:5,
                         name = c("shape_id", "shape_pt_lat", "shape_pt_lon",
                                  "shape_pt_sequence", "shape_dist_traveled", "version"),
                         type = c("TEXT", "DOUBLE", "DOUBLE", "INTEGER",
                                  "DOUBLE", "DOUBLE"),
                         notnull = 0L, dflt_value = as.logical(NA),
                         pk = c(1L, 0L, 0L, 2L, 0L, 0L),
                         stringsAsFactors = FALSE))
}
