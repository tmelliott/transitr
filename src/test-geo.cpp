#include <testthat.h>

#include "geo.h"

context("Geographical functions") {

  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("distanceEarth computes the right thing") {
    expect_true(distanceEarth(-33.2345, 173.2145, -33.2346, 173.2341) - 1822.9753132 < 0.000001);
  }

}
