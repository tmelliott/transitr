library(RSQLite)
library(RPostgreSQL)

HOST <- Sys.getenv("DATABASE_HOST")

if (HOST != "") {
    c1 <- dbConnect(SQLite(), "fulldata.db")
    vps <- dbReadTable(c1, "vehicles")
    dbDisconnect(c1)

    vps$timestamp <- as.integer(vps$timestamp)

    c2 <- dbConnect(PostgreSQL(),
                    user = "tell029", host = HOST, dbname = "realtime")
    dbWriteTable(c2, "vehicles", vps, overwrite = TRUE)
    dbDisconnect(c2)
}

