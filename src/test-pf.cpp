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
