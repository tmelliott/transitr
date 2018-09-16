# transitr

[![Travis build status](https://travis-ci.org/tmelliott/transitr.svg?branch=master)](https://travis-ci.org/tmelliott/transitr)
[![codecov](https://codecov.io/gh/tmelliott/transitr/branch/master/graph/badge.svg)](https://codecov.io/gh/tmelliott/transitr)

The goals of `transitr` are to make it easy to __load GTFS data__ into a database,
construct a __transit network__ of roads and intersections,
and __model vehicles in real-time__ from an API feed to update the network
and __generate ETAs__.


# install

`transitr` is not (yet) on CRAN, so for you would need to use `devtools`:
```r
devtools::install_github('tmelliott/transitr')
```


# usage

__Still under development!__
Please wait until the `master` branch gets pushed with the first
usable release before trying to use this package.
This here is just for demonstration of what it could be like at some point
in the future.

```r
library(transitr)
library(magrittr)

## Create a database, construct network, and connect to a realtime feed
dbname <- "realtime.db"
nw <- create_gtfs("https://cdn01.at.govt.nz/data/gtfs.zip", db = dbname) %>%
    construct() %>%
    realtime_feed("https://api.at.govt.nz/v2/public/realtime/vehiclelocations",
                  with_headers("Ocp-Apim-Subscription-Key" = "mykey"),
                  response = "protobuf")

## Set the parameters and then run the model
nw %>% 
    set_parameters(n_core = 2, 
                   n_particles = 500, 
                   gps_error = 5) %>%
    model()
```

Once running, you can launch a new R session and view the shiny app:
```r
transitr::view_realtime("realtime.db")
```


# mock data server

In order to facilitate model development and checking, there's also a mock data server
in the `simulations` directory.

To install:
```bash
cd simulations
yarn 

## or if you don't use yarn
npm install
```

To start the server, you need first an archive of vehicle position feeds,
```bash
ls archive | grep vehicle | head -n 5
# vehicle_locations_20180911050001.pb
# vehicle_locations_20180911050031.pb
# vehicle_locations_20180911050102.pb
# vehicle_locations_20180911050132.pb
# vehicle_locations_20180911050201.pb

yarn start
# yarn run v1.9.4
# $ node mock_server.js
# Mock GTFS server running on port 3000!
```

Now you can run the model with the local server, which will automatically serve 
the next file with each request.
```r
## assumeing you've constructed with simulation flag:
## $ make FLAGS="-DSIMULATION"
## simulation history will be saved in a `history` directory
dir.create("history")

## set some process ID for the server to recognise (allows running multiple simulations simultaneously)
pid <- "test1"
nw <- load_gtfs("fulldata.db") %>%
    realtime_feed(sprintf("http://localhost:3000/%s/vehicle_positions", pid),
                  response = "protobuf") %>%
    set_parameters(n_core = 1,
                   n_particles = 2000,
                   gps_error = 10)

nw %>% model()
```
