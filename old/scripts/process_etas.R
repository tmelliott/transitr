library(parallel)
sim <- commandArgs(trailingOnly = TRUE)[1]
files <- list.files(file.path("simulations", sim, "etas"), pattern = "*.pb")
times <- as.integer(gsub("etas_|\\.pb", "", files))

c <- makeCluster(3L)
q <- clusterEvalQ(c, source("scripts/common.R"))
rm(q)
clusterExport(c, "sim")
x <- pbapply::pblapply(times, function(t) loadsim(sim, t), cl = c)
rm(x)
stopCluster(c)


source("scripts/common.R")
sim <- "sim001"

for (i in seq_along(times)) {
    cat("\r", i, "of", length(times))
    loadsim(sim, times[i])
}
