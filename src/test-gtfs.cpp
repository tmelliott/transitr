#include <testthat.h>
#include <string>
#include <unistd.h>
#define GetCurrentDir getcwd

#include "gtfs.h"

context("GTFS classes") {
    // run tests from tests/testthat/ directory
    char buff[FILENAME_MAX];
    GetCurrentDir( buff, FILENAME_MAX );
    std::string current_working_dir(buff);
    std::cout << "\nCurrent path is " << current_working_dir << "\n";
    std::string dbname ("auckland_gtfs.db");
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
        std::string a0n ("NZB");
        Gtfs::Agency* a0 = gtfs.find_agency (a0n);
        expect_true (a0->agency_name () == "New Zealand Bus");
    }

}
