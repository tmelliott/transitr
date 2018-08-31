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




con <- dbConnect(SQLite(), db)
vmax <- con %>% tbl("routes") %>% 
    summarize(v = max(version, na.rm = TRUE)) %>%
    collect %>% pluck("v")

shapes <- con %>% tbl("trips") %>% 
  inner_join(con %>% tbl("routes"), by = "route_id") %>% 
  inner_join(con %>% tbl("shapes"), by = "shape_id") %>%
  filter(route_type == 3 & version == vmax) %>%
  select(shape_id, shape_pt_lon, shape_pt_lat, shape_pt_sequence) %>%
  collect

stops <- con %>% tbl("stops") %>%
  filter(version == vmax) %>%
  select(stop_id, stop_lon, stop_lat) %>%
  collect

dbDisconnect(con)

shapessf <- do.call(st_sfc, tapply(1:nrow(shapes), shapes$shape_id, function(i) {
    shapes[i,] %>% arrange(shape_pt_sequence) %>% 
      select("shape_pt_lon", "shape_pt_lat") %>% 
      as.matrix %>% st_linestring
}))

plot(shapessf)
st_crs(shapessf) <- 4326

CEN <- shapes %>% select(shape_pt_lat, shape_pt_lon) %>% as.matrix %>% colMeans %>% as.numeric

shapest <- st_transform(shapessf, sprintf("+proj=eqc +lat_0=%s +lon_0=%s", CEN[1], CEN[2]))

plot(shapest)

buf <- st_buffer(shapest, 20)

cen <- c(0, 0)
wd <- 5000
plot(buf, xlim = cen[1] + wd * c(-1, 1), ylim = cen[2] + wd * c(-1, 1))

segs <- st_sfc(crs = st_crs(buf))
for (i in 1:(length(buf)-1)) {
  for (j in (i+1):length(buf)) {
    cat(sprintf(" - %i - %i       \r", i, j))
    segij <- st_intersection(buf[[i]], buf[[j]])
    if (length(segij)) {
      for (k in 1:length(segij)) {
        area <- segij[[1]] %>% st_polygon %>% st_area
        area / 40
      }
      break()
      #segs <- c(segs, list(segij))
    }
  }
}

segs <- do.call(st_sfc, segs)


plot(do.call(st_sfc, segs), xlim = cen[1] + wd * c(-1, 1), ylim = cen[2] + wd * c(-1, 1))


### 
sht <- st_transform(shapessf, sprintf("+proj=eqc +lat_0=%s +lon_0=%s", CEN[1], CEN[2]))
stps <- stops[,2:3] %>% as.matrix %>% st_multipoint %>%
  st_sfc(crs = 4326) %>%
  st_transform(crs = st_crs(sht))


s1 <- st_sfc(sht[[1]], crs = st_crs(sht))
plot(s1)
plot(stps, add= T)



s2 <- st_buffer(st_sfc(sht[[6]], crs = st_crs(sht)), 20)

xx <- st_intersection(s1, s2)

apply(as.matrix(s1), 1, function(x) st_distance(st_sfc(st_point(x), crs = 4326), s2))



