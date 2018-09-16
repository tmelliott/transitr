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
