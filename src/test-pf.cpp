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

    test_that("Segment index is returned correctly") {
        std::vector<Gtfs::ShapeSegment> segments;
        segments.emplace_back ();
        segments.back ().distance = 0;
        segments.emplace_back ();
        segments.back ().distance = 1000;
        segments.emplace_back ();
        segments.back ().distance = 2000;
        segments.emplace_back ();
        segments.back ().distance = 3000;

        expect_true (Gtfs::find_segment_index (0, &segments) == 0);
        expect_true (Gtfs::find_segment_index (100, &segments) == 0);
        expect_true (Gtfs::find_segment_index (1000, &segments) == 1);
        expect_true (Gtfs::find_segment_index (1100, &segments) == 1);
        expect_true (Gtfs::find_segment_index (3000, &segments) == 3);
        expect_true (Gtfs::find_segment_index (4000, &segments) == 3);
    }

    test_that("Segment index is returned correctly even if circular goes to 0") {
        std::vector<Gtfs::ShapeSegment> seg2;
        seg2.emplace_back ();
        seg2.back ().distance = 0;
        seg2.emplace_back ();
        seg2.back ().distance = 1000;
        seg2.emplace_back ();
        seg2.back ().distance = 2000;
        seg2.emplace_back ();
        seg2.back ().distance = 0;

        expect_true (Gtfs::find_segment_index (0, &seg2) == 0);
        expect_true (Gtfs::find_segment_index (100, &seg2) == 0);
        expect_true (Gtfs::find_segment_index (1000, &seg2) == 0);
        expect_true (Gtfs::find_segment_index (1100, &seg2) == 0);
        expect_true (Gtfs::find_segment_index (3000, &seg2) == 0);
    }
}

context ("Vehicle states") {
    std::string dbname ("auckland_gtfs.db");
    Gtfs::Gtfs gtfs (dbname);
    Gtfs::par par;

    RNG rng (10);

    std::string vid ("test");
    Gtfs::Vehicle v (vid, &par);

    std::string t ("1141160875-20190613111133_v80.31");
    Gtfs::Trip* t0 = gtfs.find_trip (t);

    test_that ("Vehicle loads OK with a trip") {
        expect_true (v.trip () == nullptr);
        v.set_trip (t0);
        expect_false (v.trip () == nullptr);
    }

    test_that ("Vehicle initializes to zero with position update at stop 1") {
        uint64_t ts = 1562034647;
        Gtfs::Stop* s1 = t0->stops ().at (0).stop;
        v.add_event (Gtfs::Event (ts, Gtfs::EventType::gps, t, s1->stop_position ()));
        v.update (&gtfs);
        expect_true (v.get_events ().size () == 1);

        v.mutate (rng, &gtfs);
        expect_true (v.state ()->size () == par.n_particles);
    }
}
