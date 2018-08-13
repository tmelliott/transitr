library(tidyverse)

vid <- "7D96"

ph <- read_csv(sprintf("vehicle_%s.csv", vid),
               col_names = c("timestamp", "distance", "speed", "lat", "lon")) %>%
    mutate(timestamp = as.POSIXct(timestamp, origin = "1970-01-01")) 
gridExtra::grid.arrange(
    ggplot(ph, aes(lon, lat)) + geom_point(),
    ggplot(ph, aes(timestamp, distance)) + geom_point(),
    heights = c(3, 1))
