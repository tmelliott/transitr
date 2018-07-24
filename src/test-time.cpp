#include <testthat.h>

#include "time.h"

context("Time functions") {

  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  test_that("Time values are correct") {
    Time tx = Time (10, 20, 30); // 10:20:30
    expect_true(tx.hour () == 10);
    expect_true(tx.minute () == 20);
    expect_true(tx.second () == 30);
  }

}
