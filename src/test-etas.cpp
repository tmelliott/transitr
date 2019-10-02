#include <testthat.h>
#include <iostream>
#include <string>

#include "gtfs.h"

context("Arrival time estimation") {
    // first need GTFS object
    std::string dbname ("auckland_gtfs.db");
    Gtfs::Gtfs gtfs (dbname);
    Gtfs::par par;
    par.nw_measurement_error = 5;
    RNG rng (10);

    test_that("Trip states initialized correctly") {
        // create a trip and initialize it .. 
        std::string t ("1141160875-20190613111133_v80.31");
        Gtfs::Trip* t0 = gtfs.find_trip (t);

        uint64_t ts = 1562034647;

        t0->update (ts, rng);

        expect_true (0);
    }
}
