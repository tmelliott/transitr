# transitr

[![Travis build status](https://travis-ci.org/tmelliott/transitr.svg?branch=develop)](https://travis-ci.org/tmelliott/transitr)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/tmelliott/transitr?branch=develop&svg=true)](https://ci.appveyor.com/project/tmelliott/transitr)
[![codecov](https://codecov.io/gh/tmelliott/transitr/branch/develop/graph/badge.svg)](https://codecov.io/gh/tmelliott/transitr)

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
nw <- create_nw("https://cdn01.at.govt.nz/data/gtfs.zip") %>%
   construct() %>%
   connect("https://api.at.govt.nz/v2/public/realtime/vehiclelocations") %>%
   with_headers("Ocp-Apim-Subscription-Key" = "mykey")

## Set the model going - hopefully this will spawn a new child R process
## that continually runs in the background ... (until you kill it)
## note: n.particles should be bigger than 500,
##       and cores should be as many as you have spare
nw %>% model(cores = 2, n.particles = 500)

## Once running, you can view ETAs for buses arriving at a stop
nw %>% stop('1529') %>% etas()

## or, maybe at some point this will open up a shiny app to explore all the things
## in real time???
nw %>% view()

## and when finished
nw %>% close()
## will close the child process and free up your computers resources
```