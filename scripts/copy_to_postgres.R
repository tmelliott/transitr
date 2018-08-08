library(RSQLite)
library(RPostgreSQL)

HOST <- Sys.getenv("DATABASE_HOST")

if (HOST != "") {
    c1 <- dbConnect(SQLite(), "fulldata.db")
    vps <- dbReadTable(c1, "vehicles")
    dbDisconnect(c1)

    vps$timestamp <- as.POSIXct(as.integer(vps$timestamp), origin = "1970-01-01")

    c2 <- dbConnect(PostgreSQL(),
                    user = "tell029", host = HOST, dbname = "realtime")
    dbGetQuery(c2, paste(
        "INSERT INTO vehicles_history",
        "  SELECT * FROM vehicles"
        "ON CONFLICT DO NOTHING"))
    dbGetQuery(c2, "DELETE FROM vehicles")
    dbWriteTable(c2, "vehicles", vps, append = TRUE, row.names = FALSE)
    dbDisconnect(c2)
}

