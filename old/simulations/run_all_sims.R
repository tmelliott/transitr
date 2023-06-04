simdirs <- list.files(pattern = "sim_")

for (simdir in simdirs) {
    if (file.exists(file.path(simdir, "timings.csv"))) next()
    # system(sprintf("cd ~/transitr && make simulation SIM=%s", simdir))
    system(sprintf("cd .. && make simulation SIM=%s", simdir))
}

