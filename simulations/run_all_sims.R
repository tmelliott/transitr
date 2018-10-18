simdirs <- list.files(pattern = "sim_")

for (simdir in simdirs)
    system(sprintf("cd ~/transitr && make simulation SIM=%s", simdir))

