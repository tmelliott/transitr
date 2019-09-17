#include <testthat.h>
#include <iostream>
#include <string>

#include "gtfs.h"


context("Network filter") {
    std::string dbname ("auckland_gtfs.db");
    Gtfs::Gtfs gtfs (dbname);
    Gtfs::par par;
    RNG rng (10);

    test_that ("Segment states initialize") {
        Gtfs::Segment* seg = gtfs.find_segment (1);
        expect_true (seg->segment_id () == 1);
        // the segment initialisation is done it a call to update ()
        std::cout << "\n----\nBeginning segment test 1\n";
        // seg->update (&par, &gtfs);
    //     expect_true (seg->is_loaded ());
        std::cout << "\n----\nComplete\n\n";
    }
}
