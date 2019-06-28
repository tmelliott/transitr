#include <testthat.h>
#include <string>

#include "gtfs.h"

context("GTFS classes") {
    // run tests from tests/testthat/ directory
    std::string dbname ("auckland_gtfs.db");
    Gtfs::Gtfs gtfs (dbname);

    test_that("database connects and loads templates") {
        expect_true (gtfs.agencies ().size () == 16);
        expect_true (gtfs.routes ().size () == 2079);
        expect_true (gtfs.trips ().size () == 15);
        expect_true (gtfs.shapes ().size () == 15);
        expect_true (gtfs.stops ().size () == 356);
        expect_true (gtfs.calendar ().size () == 15);
        expect_true (gtfs.nodes ().size () == 356);
    }

    test_that("objects load on request") {
        std::string t ("1141160875-20190613111133_v80.31");
        Gtfs::Trip* t0 = gtfs.find_trip (t);
        expect_true (t0->trip_id () == t);
        expect_true (t0->trip_headsign () == "Britomart");
        expect_true (t0->route ()->route_short_name () == "70");
        expect_true (t0->route ()->agency ()->agency_name () == "Howick and Eastern");


        expect_true (t0->shape ()->path ().size () == 1019);
        expect_true (t0->shape ()->nodes ().size () == 47);
        // expect_true (t0->shape ()->segments ().size () == 1);
        expect_true (t0->calendar ()->monday ());
        // expect_true (t0->calendar ()->exceptions ().size () == 0);
        expect_true (t0->stops ().size () == 47);
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
        std::vector<double> ds {0.0, 1000.0, 2000.0, 10000.0, 100000.0};
        double dx;
        latlng p;
        double slen = s->path ().back ().distance;
        std::cout << "\n shape_id = " << s->shape_id () << ", length = "
            << slen << "\n";
        for (auto d : ds)
        {
            p = s->coordinates_of (d);
            dx = s->distance_of (p);
            std::cout << "\n " << d << " -> " << dx;
            expect_true (fabs (fmin (d, slen) - dx) < 1);
        }
        std::cout << "\n";
    }

    test_that("Nodes load") {
        Gtfs::Node* n = gtfs.find_node (1);

    }

}
