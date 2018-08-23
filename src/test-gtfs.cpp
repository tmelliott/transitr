#include <testthat.h>
#include <string>

#include "gtfs.h"

context("GTFS classes") {
    // run tests from tests/testthat/ directory
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
        std::string t ("1141101952-20180702170310_v67.28");
        Gtfs::Trip* t0 = gtfs.find_trip (t);
        expect_true (t0->trip_id () == t);
        expect_true (t0->trip_headsign () == "Britomart");
        expect_true (t0->route ()->route_short_name () == "70");
        expect_true (t0->route ()->agency ()->agency_name () == "Howick and Eastern");


        expect_true (t0->shape ()->path ().size () == 1030);
        // expect_true (t0->shape ()->segments ().size () == 1);
        expect_true (t0->calendar ()->monday ());
        // expect_true (t0->calendar ()->exceptions ().size () == 0);
        expect_true (t0->stops ().size () == 48);
    }

    test_that("Invalid requests return nullptrs") {
        std::string ne ("nonexistent");
        expect_true (gtfs.find_trip (ne) == nullptr);
        expect_true (gtfs.find_route (ne) == nullptr);
        expect_true (gtfs.find_agency (ne) == nullptr);
        expect_true (gtfs.find_shape (ne) == nullptr);
        expect_true (gtfs.find_calendar (ne) == nullptr);
    }

    test_that("Shape distance functions return the correct values") {
        Gtfs::Shape* s = &(gtfs.shapes ().begin (0)->second);
        double d (1000);
        latlng p = s->coordinates_of (d);
        expect_true (fabs (s->distance_of (p) - d) < 0.001);
    }

}
