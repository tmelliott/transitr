#include <testthat.h>
#include <string>

#include "gtfs.h"

context("Particle filter helper functions") {
    test_that("Stop index is returned correctly") {
        std::vector<Gtfs::StopTime> stops;
        stops.emplace_back ();
        stops.back ().distance = 0;
        stops.emplace_back ();
        stops.back ().distance = 1000;
        stops.emplace_back ();
        stops.back ().distance = 2000;
        stops.emplace_back ();
        stops.back ().distance = 3000;

        expect_true (Gtfs::find_stop_index (0, &stops) == 0);
        expect_true (Gtfs::find_stop_index (100, &stops) == 0);
        expect_true (Gtfs::find_stop_index (1000, &stops) == 1);
        expect_true (Gtfs::find_stop_index (1100, &stops) == 1);
        expect_true (Gtfs::find_stop_index (3000, &stops) == 3);
        expect_true (Gtfs::find_stop_index (4000, &stops) == 3);
    }
}
