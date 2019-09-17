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
            expect_true (exp (p->get_ll ()) > 0);
        }
    }

    d = stops.at (2).distance - stops.at (0).distance;
    uint64_t ts4 = ts3 + round (d / 14.0);
    v.add_event (Gtfs::Event (ts4, Gtfs::EventType::departure, t, 2));
    v.update (&gtfs);
    v.mutate (rng, &gtfs);

    // std::cout << "\nSeg travel time (0->1): "
    //     << v.segment_travel_time (0) << "\n";

    test_that("Travel times estimated") {
        expect_true (v.segment_travel_time (0) > 0);
        expect_true (v.state ()->size () == 10);
        for (auto p = v.state ()->begin (); p != v.state ()->end (); ++p)
        {
            expect_true (p->get_departure_time (2) > ts);
        }
    }


    int M = stops.size () - 1;
    d = stops.at (M).distance - stops.at (2).distance;
    uint64_t ts5 = ts4 + round (d * 2) + 100 * M;
    v.add_event (Gtfs::Event (ts5, Gtfs::EventType::arrival, t, M));
    v.update (&gtfs);
    v.mutate (rng, &gtfs);

    test_that ("All travel times get computed ...") {
        for (int i=0; i<M-1; i++)
        {
            expect_true (v.segment_travel_time (i) > 0);
        }
    }
}


context("Vehicle mutate/update from GPS obs") {
    std::string dbname ("auckland_gtfs.db");
    Gtfs::Gtfs gtfs (dbname);
    Gtfs::par par;
    par.n_particles = 1000;
    par.gps_error = 20.0;
    par.system_noise = 2;

    RNG rng (10);

    std::string vid ("test");
    Gtfs::Vehicle v (vid, &par);

    std::string t ("1141160875-20190613111133_v80.31");
    Gtfs::Trip* t0 = gtfs.find_trip (t);
    Gtfs::Shape* shape = t0->shape ();
    double Dmax = shape->path ().back ().distance;
    uint64_t ts = 1562034647;
    
    std::vector<double> x {0.1, 0.18, 0.2, 0.3, 0.35, 0.5, 0.8, 0.9, 0.95};
    // std::vector<double> xdot {12,   14,  10,  10,   20,  18,  17,  10};
    double xdot = 12;
    // let's assume the bus never stops ...
    double xi = Dmax * x.at (0);
    latlng pos = shape->coordinates_of (xi);

    v.add_event (Gtfs::Event (ts, Gtfs::EventType::gps, t, pos));
    v.update (&gtfs);
    v.mutate (rng, &gtfs);

    double xi1;
    int delta;
    test_that ("Vehicle initialized OK from GPS obs") {
        expect_true (v.state ()->size () == par.n_particles);

        for (int i=1; i<x.size (); i++)
        {
            xi1 = xi;
            xi = Dmax * x.at (i);
            delta = (xi - xi1) / xdot; //.at (i-1);
            pos = shape->coordinates_of (xi);
            ts += delta;
            v.add_event (Gtfs::Event (ts, Gtfs::EventType::gps, t, pos));
            v.update (&gtfs);
            v.mutate (rng, &gtfs);
        }
        // finally, STU for last stop
        delta = (Dmax - xi) / xdot;
        ts += delta;
        v.add_event (
            Gtfs::Event (
                ts, 
                Gtfs::EventType::arrival, 
                t, 
                t0->stops ().size () - 1
            )
        );
        v.update (&gtfs);
        v.mutate (rng, &gtfs);
        for (int i=5; i<t0->shape ()->segments ().size () - 1; i++)
        {
            std::cout << "\n seg " << i << " = " << v.segment_travel_time (i);
            expect_true (v.segment_travel_time (i) > 0);
        }
    }


}
