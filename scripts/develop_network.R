library(tidyverse)
library(RSQLite)
library(dbplyr)

db <- "fulldata.db"
con <- dbConnect(SQLite(), db)

vmax <- con %>% tbl("routes") %>% 
    summarize(v = max(version, na.rm = TRUE)) %>%
    collect %>% pluck("v")
routes <- con %>% tbl("routes") %>%
    filter(route_short_name %in% c("982", "983") &
           route_long_name %like% "Gulf Harbour To HC Station%" &
           version == vmax) %>%
    inner_join(con %>% tbl("trips"), by = "route_id") %>%
    select(route_id, route_short_name, route_long_name, shape_id) %>%
    distinct %>% collect

shapes <- con %>% tbl("shapes") %>%
    filter(shape_id %in% routes$shape_id) %>%
    arrange(shape_id, shape_pt_sequence) %>%
    select(shape_id, shape_pt_lat, shape_pt_lon, shape_pt_sequence) %>%
    collect

latr <- extendrange(shapes$shape_pt_lat)
lonr <- extendrange(shapes$shape_pt_lon)
stops <- con %>% tbl("stops") %>%
    filter(stop_lat > latr[1] & stop_lat < latr[2] &
           stop_lon > lonr[1] & stop_lon < lonr[2] &
           version == vmax) %>%
    select(stop_lon, stop_lat) %>%
    collect

dbDisconnect(con)

ggplot(shapes, aes(shape_pt_lon, shape_pt_lat, 
                   group = shape_id, colour = shape_id)) +
    geom_path()




## should end up with 15 segments; 12 nodes
path1 <- shapes %>% 
    filter(shape_id == routes$shape_id[routes$route_short_name == "983"]) %>%
    arrange(shape_pt_sequence) %>% 
    delete_stops(stops)
path2 <- shapes %>% 
    filter(shape_id == routes$shape_id[routes$route_short_name == "982"]) %>%
    arrange(shape_pt_sequence) %>% 
    delete_stops(stops)

existing <- path1[, 3:2]
new <- path2[, 3:2]

# ds <- merge_paths(new, existing)
# ggplot(NULL, aes(seq_along(ds), ds, color = ds < 5)) + geom_point()

# p <- ggplot(path1, aes(shape_pt_lon, shape_pt_lat)) +
#     geom_path() +
#     geom_point(aes(colour = d < 5), data = path2 %>% mutate(d = ds) %>% arrange(d))
# p
# # p + xlim(174.73, 174.74) + ylim(-36.63, -36.625)

# di <- ds < 5
# d1 <- d2 <- integer(nrow(path2))
# for (i in 2:length(d1)) if (di[i] == di[i-1]) d1[i] <- d1[i-1] + 1
# for (i in (length(d2)-1):1) if (di[i] == di[i+1]) d2[i] <- d2[i+1] + 1

# n.thres <- 5
# if (any(d1 + d2 < n.thres)) {
#     di[d1 + d2 < n.thres] <- !di[d1 + d2 < n.thres]
#     d1 <- d2 <- integer(nrow(path2))
#     for (i in 2:length(d1)) if (di[i] == di[i-1]) d1[i] <- d1[i-1] + 1
#     for (i in (length(d2)-1):1) if (di[i] == di[i+1]) d2[i] <- d2[i+1] + 1    
# }


di <- merge_paths(new, existing)

ggplot(path1 %>% mutate(d = existing.seg) %>% filter(!is.na(d)), 
       aes(shape_pt_lon, shape_pt_lat, color = as.factor(d))) +
    geom_path() +
    geom_path(aes(colour = as.factor(d)), data = path2 %>% mutate(d = new.seg + max(existing.seg, na.rm = TRUE)))


ggplot(segs, aes(shape_pt_lon, shape_pt_lat, colour = as.factor(i))) +
    geom_path()
