paste_nl <- function(...) paste(..., sep = "\n")

db_connect <- function(x) UseMethod('db_connect')
db_connect.default <- function(x) RSQLite::dbConnect(RSQLite::SQLite(), x)
db_connect.trgtfs <- function(x) db_connect(x$database)
db_close <- function(con) RSQLite::dbDisconnect(con)
