library(tidyverse)
library(RProtoBuf)
library(RSQLite)
library(dbplyr)

load("tripupdates.rda")
db <- "fulldata.db"

## Filter out the westgate buses: 125X and 111
con <- dbConnect(SQLite(), db)
routes <- con %>% tbl("routes") %>%
    filter(route_short_name %in% c('125X', '111')) %>%
    filter(route_long_name %like% 'City Centre%' | route_long_name == 'Royal Heights Loop') %>%
    collect() %>%
    group_by(route_short_name) %>%
    summarize(route_id = first(route_id)) %>%
    mutate(route_id = gsub("-.*", "", route_id))

tripupdates <- tripupdates %>%
    filter(grepl(paste0('^', routes$route_id, collapse = "|"), route_id))

arr <- tripupdates %>%
    filter(grepl(paste0('^', routes$route_id[2]), route_id) & stop_sequence == 2)

dep <- tripupdates %>%
    filter(grepl(paste0('^', routes$route_id[1]), route_id) & stop_sequence == 1)

# for each arr, find next dep
