#include <testthat.h>
#include <iostream>
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
    par.n_particles = 10;

    RNG rng (10);

    std::string vid ("test");
    Gtfs::Vehicle v (vid, &par);

    std::string t ("1141160875-20190613111133_v80.31");
    Gtfs::Trip* t0 = gtfs.find_trip (t);
    uint64_t ts = 1562034647;
    Gtfs::Stop* s1 = t0->stops ().at (0).stop;

    test_that ("Vehicle loads OK with a trip") {
        expect_true (v.trip () == nullptr);
        v.set_trip (t0);
        expect_false (v.trip () == nullptr);
    }

    test_that ("Vehicle initializes OK with position update at stop 1") {
        // this timestamp holds no significance, other than it was the time
        // when I needed a timestamp!
        v.add_event (Gtfs::Event (ts, Gtfs::EventType::gps, t, s1->stop_position ()));
        v.update (&gtfs);
        expect_true (v.get_events ().size () == 1);

        v.mutate (rng, &gtfs);
        expect_true (v.state ()->size () == par.n_particles);
        for (auto p = v.state ()->begin (); p != v.state ()->end (); ++p)
        {
            expect_true (p->get_distance () <= 100);
            expect_true (
                p->get_travel_times ().size () == 
                    t0->shape ()->segments ().size ()
            );
        }
    }

    std::vector<Gtfs::StopTime>& stops = t0->stops ();
    test_that("Arrival and departure observations at stop 2 get dealt with correctly") {
        // another timestamp, I suppose
        // we pretend that the previous observation was at STOP 1
        // and bus travels at 15m/s to next stop, and waits for 20 seconds
        double d = stops.at (1).distance - stops.at (0).distance;
        uint64_t ts2, ts3;
        ts2 = ts + round(d / 15.0);
        int dwell = 20;
        ts3 = ts2 + dwell;
        v.add_event (Gtfs::Event (ts2, Gtfs::EventType::arrival, t, 1));
        v.add_event (Gtfs::Event (ts3, Gtfs::EventType::departure, t, 1));
        v.update (&gtfs);
        expect_true (v.get_events ().size () == 2);

        // now mutate that 
        v.mutate (rng, &gtfs);
        for (auto p = v.state ()->begin (); p != v.state ()->end (); ++p)
        {
            // expect_true(p->get_arrival_time (1) > 0);
        }
    }
}

context("Particle functions") {
    std::string dbname ("auckland_gtfs.db");
    Gtfs::Gtfs gtfs (dbname);
    Gtfs::par par;
    par.n_particles = 10;

    RNG rng (10);

    std::string vid ("test");
    Gtfs::Vehicle v (vid, &par);

    std::string t ("1141160875-20190613111133_v80.31");
    Gtfs::Trip* t0 = gtfs.find_trip (t);
    uint64_t ts = 1562034647;
    std::vector<Gtfs::StopTime> stops = t0->stops ();

    v.add_event (Gtfs::Event (ts, Gtfs::EventType::arrival, t, 0));
    v.update (&gtfs);
    v.mutate (rng, &gtfs);

    test_that("Particles initialized at stop with arrival update") {
        expect_true (v.state ()->size () == 10);
        for (auto p = v.state ()->begin (); p != v.state ()->end (); ++p)
        {
            expect_true (p->get_distance () == 0);
            expect_true (p->get_stop_index () == 0);
        }
    }

    double d = stops.at (1).distance - stops.at (0).distance;
    uint64_t ts2, ts3;
    ts2 = ts + round(d / 12.0);
    int dwell = 20;
    ts3 = ts2 + dwell;

    test_that("Particle travel function behaves") {
        Gtfs::Event e2 (ts2, Gtfs::EventType::arrival, t, 1);
        Gtfs::Event e3 (ts3, Gtfs::EventType::departure, t, 1);
        expect_true (v.state ()->size () == 10);
        v.override_timestamp (ts2);
        for (auto p = v.state ()->begin (); p != v.state ()->end (); ++p)
        {
            expect_true (p->get_distance () == 0);
            p->travel (ts2 - ts, e2, rng);

            expect_true (p->get_distance () >= 0);
            expect_true (p->get_distance () < 30 * (ts2 - ts));

            // travel time is extrapolated forward, even if it doesn't reach stop
            expect_true (p->get_arrival_time (1) > ts);
        }

        v.override_timestamp (ts3);
        for (auto p = v.state ()->begin (); p != v.state ()->end (); ++p)
        {
            p->travel (ts3 - ts2, e3, rng);
            expect_true (p->get_departure_time (1) > ts);
        }
    }
}

context("Vehicle mutate/update") {
    std::string dbname ("auckland_gtfs.db");
    Gtfs::Gtfs gtfs (dbname);
    Gtfs::par par;
    par.n_particles = 10;

    RNG rng (10);

    std::string vid ("test");
    Gtfs::Vehicle v (vid, &par);

    std::string t ("1141160875-20190613111133_v80.31");
    Gtfs::Trip* t0 = gtfs.find_trip (t);
    uint64_t ts = 1562034647;
    std::vector<Gtfs::StopTime> stops = t0->stops ();

    v.add_event (Gtfs::Event (ts, Gtfs::EventType::arrival, t, 0));
    v.update (&gtfs);
    v.mutate (rng, &gtfs);

    test_that("Particle weights start off OK") {
        expect_true (v.state ()->size () == 10);
        for (auto p = v.state ()->begin (); p != v.state ()->end (); ++p)
        {
            expect_true (p->get_weight () == 0.1);
        }
    }

    double d = stops.at (1).distance - stops.at (0).distance;
    uint64_t ts2, ts3;
    ts2 = ts + round(d / 12.0);
    v.add_event (Gtfs::Event (ts2, Gtfs::EventType::arrival, t, 1));
    v.update (&gtfs);
    v.mutate (rng, &gtfs);

    test_that("Particles get resampled?") {
        expect_true (v.state ()->size () == 10);
        for (auto p = v.state ()->begin (); p != v.state ()->end (); ++p)
        {
            std::cout << "\n l(y|p) = " << exp (p->get_ll ());
            // expect_true (exp (p->get_ll ()) > 0);
        }
    }

    int dwell = 20;
    ts3 = ts2 + dwell;
    v.add_event (Gtfs::Event (ts3, Gtfs::EventType::departure, t, 1));
    v.update (&gtfs);
    v.mutate (rng, &gtfs);

    test_that("Particles get resampled (2)?") {
        expect_true (v.state ()->size () == 10);
        for (auto p = v.state ()->begin (); p != v.state ()->end (); ++p)
        {
            std::cout << "\n l(y|p) = " << exp (p->get_ll ());
            expect_true (exp (p->get_ll ()) > 0);
        }
        // expect_true (1 == 2);
    }
}
