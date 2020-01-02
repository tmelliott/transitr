library(tidyverse)

eta_quantiles <-
    read_csv("simulations/sim002/eta_quantiles.csv",
        col_names = c("trip_id", "vehicle_id", "stop_sequence", "timestamp", "eta", "quantile"),
        col_types = "cciiin"
    )


i <- 1800000
T <- eta_quantiles$trip_id[i]
S <- eta_quantiles$stop_sequence[i]
TS <- eta_quantiles$timestamp[i]

eta1 <- eta_quantiles %>%
    filter(trip_id == T & stop_sequence == S & timestamp == TS)

ggplot(eta1, aes(eta, quantile)) +
    geom_path()
