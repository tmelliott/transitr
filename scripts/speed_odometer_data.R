## Step 1: load historical data set
SIM_DATE <- "2020-05-27"
if (!file.exists(sprintf("simulations/archive_%s.zip", gsub("-", "_", SIM_DATE)))) {
    system(
        sprintf(
            "scp tom@%s:/mnt/storage/history/%s/archive_%s.zip simulations/archive_%s.zip",
            Sys.getenv("pi_ip"),
            gsub("-", "/", SIM_DATE),
            gsub("-", "_", SIM_DATE),
            gsub("-", "_", SIM_DATE)
        )
    )
}
unlink("simulations/archive", TRUE, TRUE)
unzip(
    sprintf("simulations/archive_%s.zip", gsub("-", "_", SIM_DATE)),
    exdir = "simulations/archive"
)

## Step 2: ensure database populated


## Step 3: launch mock-server
system("cd simulations && nohup yarn start &")

## Step 4: run program


## Step 5: analyse results
