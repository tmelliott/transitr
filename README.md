
<!-- README.md is generated from README.Rmd. Please edit that file -->

# transitr

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/tmelliott/transitr/branch/develop/graph/badge.svg)](https://app.codecov.io/gh/tmelliott/transitr?branch=develop)
[![R-CMD-check](https://github.com/tmelliott/transitr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tmelliott/transitr/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

The goals of `transitr` are to make it easy to **load GTFS data** into a
database, construct a **transit network** of roads and intersections,
and **model vehicles in real-time** from an API feed to update the
network and **generate ETAs**.

## Installation

You can install the development version of transitr from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tmelliott/transitr")
```

## Example

**Still under development\!** This here is just for demonstration of
what it could be like at some point in the future.

``` r
library(transitr)

## Create a database, construct network, and connect to a realtime feed
dbname <- "realtime.db"
nw <- create_gtfs("https://cdn01.at.govt.nz/data/gtfs.zip", db = dbname) |>
    construct() |>
    realtime_feed("https://api.at.govt.nz/v2/public/realtime/vehiclelocations",
                  with_headers("Ocp-Apim-Subscription-Key" = "mykey"),
                  response = "protobuf")

## Set the parameters and then run the model
nw |>
    set_parameters(n_core = 2,
                   n_particles = 500,
                   gps_error = 5) |>
    model()
```
