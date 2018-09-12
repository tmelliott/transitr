library(tidyverse)

whatOrder <- c("number of vehicles", "loading vehicle positions", 
               "updating vehicle information", "updating vehicle states", 
               "predicting ETAs", "writing ETAs to protobuf feed")
view <- function(sim) {
    timings <- read_csv(file.path("simulations", sim, "timings.csv"))
    timings <- timings %>% bind_rows(timings %>% group_by(iteration) %>% 
            summarize(timestamp = mean(timestamp), nvehicles = mean(nvehicles)) %>%
            mutate(what = "number of vehicles",
                   cpu = nvehicles, wall = nvehicles)
            ) %>%
        mutate(what = fct_relevel(what, whatOrder),
               timestamp = as.POSIXct(timestamp, origin = "1970-01-01"))

    ggplot(timings, aes(timestamp, cpu)) +
        geom_line() +
        facet_grid(what~., scales = "free_y")
}

view("sim000")
