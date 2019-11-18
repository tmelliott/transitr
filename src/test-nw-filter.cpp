#include <testthat.h>
#include <iostream>
#include <string>

#include "gtfs.h"


context("Network filter") {
    std::string dbname ("auckland_gtfs.db");
    Gtfs::Gtfs gtfs (dbname);
    Gtfs::par par;
    par.nw_measurement_error = 5;
    RNG rng (10);

    test_that ("Segment states initialize") {
        Gtfs::Segment* seg = gtfs.find_segment (1);
        expect_true (seg->segment_id () == 1);
        // the segment initialisation is done it a call to update ()
        std::cout << "\n----\nBeginning segment test 1\n";
        seg->update (&par, &gtfs);
        expect_true (seg->is_loaded ());

        expect_true (seg->speed () == 30. / 3.6);
        std::cout << "\n - segment length = "
            << seg->length () << "m"
            << "\n - travel time state = "
            << seg->speed ()
            << " (" << seg->uncertainty () << ") m/s";

        double obs, err;
        obs = 40;
        err = 3;
        uint64_t ts = 1562034647;

        // IKF
        double B, P, U, u, I, i;
        B = seg->speed ();
        P = 5.;

        std::cout << "\n\n - pretend observation: "
            << obs << " \u00B1 " << err << " \u00B1 " << seg->state_var ();

        expect_true (seg->get_data ().size () == 0);
        seg->push_data (obs, err, ts);
        expect_true (seg->get_data ().size () == 1);
        seg->update (ts);
        // data should be cleared
        expect_true (seg->get_data ().size () == 0);

        U = 1 / P;
        u = B / P;
        err = par.nw_measurement_error + pow(seg->state_var (), 2);
        I = 1 / err;
        i = obs / err;
        U += I;
        u += i;
        P = 1 / U;
        B = u / U;

        expect_true (seg->speed () == B);

        std::cout << "\n----\nComplete\n\n";

        // expect_true (1 == 0);
    }
}
