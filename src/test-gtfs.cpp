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
        expect_true (t0->stops ().size () == 47);
        // expect_true (t0->shape ()->segments ().size () == 1);
        expect_false (t0->calendar ()->monday ());
        expect_true (t0->calendar ()->sunday ());
        // expect_true (t0->calendar ()->exceptions ().size () == 0);
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
        for (auto d : ds)
        {
            p = s->coordinates_of (d);
            dx = s->distance_of (p);
            expect_true (fabs (fmin (d, slen) - dx) < 1);
        }
    }

    test_that("Nodes load") {
        Gtfs::Node* n = gtfs.find_node (1);
        expect_true (n->node_type () == 0);
    }

    test_that("Stop and node distances are equal") {
        std::string t ("1141160875-20190613111133_v80.31");
        Gtfs::Trip* t0 = gtfs.find_trip (t);
        std::vector<Gtfs::StopTime> stops = t0->stops ();
        std::vector<Gtfs::ShapeNode> nodes = t0->shape ()->nodes ();
        expect_true (stops.size () == nodes.size ());
        for (int i=0; i<stops.size (); i++)
        {
            expect_true (stops.at (i).distance == nodes.at (i).distance);
        }
    }

    test_that("Shape segments are loaded correctly") {
        Gtfs::Shape* s = &(gtfs.shapes ().begin (0)->second);
        std::vector<Gtfs::ShapeSegment> segs = s->segments ();
        expect_true (segs.size () == s->nodes ().size () - 1);
    }

    test_that("Path length is valid") {
        for (auto t = gtfs.trips ().begin (); t != gtfs.trips ().end (); ++t)
        {
            Gtfs::Trip* t0 = &(t->second);
            std::vector<Gtfs::StopTime> stops = t0->stops ();
            Gtfs::Shape* s = t0->shape ();
            std::cout << "\n";
            double d (0);
            latlng *p1, *p2;
            p2 = &(s->path ().at (0).pt);
            for (int i=1; i<s->path ().size (); i++)
            {
                p1 = p2;
                p2 = &(s->path ().at (i).pt);

                d += distanceEarth (*p1, *p2);
                expect_true (d == s->path ().at (i).distance);
            }

            std::cout << "\n stops len = " << stops.back ().distance;
            std::cout << "\n nodes len = " << s->nodes ().back ().distance;
            std::cout << "\n path len = " << s->path ().back ().distance;
            expect_true (s->path ().back ().distance - stops.back ().distance < 0.1);
            expect_true (s->path ().back ().distance - s->nodes ().back ().distance < 0.5);
        }
    }
}


context ("Road network") {
    std::string dbname ("auckland_gtfs.db");
    Gtfs::Gtfs gtfs (dbname);

    Gtfs::par par;

    test_that("Initialising road network works") {
        for (auto it = gtfs.segments ().begin (); it != gtfs.segments ().end (); ++it)
        {
            it->second.update (&par, &gtfs);
            expect_true (it->second.travel_time () == it->second.length () / 10.0);
        }
    }
}


