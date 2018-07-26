#include <testthat.h>
#include <string>

#include "gtfs.h"

context("GTFS classes") {
    // run tests from tests/testthat/ directory
    std::string dbname ("../../inst/extdata/auckland_gtfs.db");
    Gtfs::Gtfs gtfs (dbname);

    test_that("database connects and loads templates") {
        expect_true (gtfs.agencies ().size () == 17);
        expect_true (gtfs.routes ().size () == 2234);
        expect_true (gtfs.trips ().size () == 17);
        expect_true (gtfs.shapes ().size () == 17);
        expect_true (gtfs.stops ().size () == 369);
        expect_true (gtfs.calendar ().size () == 17);
    }

    test_that("objects load on request") {
        std::string a0n ("Go West");
        Gtfs::Agency* a0 = gtfs.find_agency (a0n);
        std::cout << a0->agency_name ();
        expect_true (a0->agency_name () == a0n);
    }

}
