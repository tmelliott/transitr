#include <testthat.h>

#include "time.h"

context("Time functions") {

  test_that("Time as three integers") {
    Time tx = Time (10, 20, 30); // 10:20:30
    expect_true (tx.hour () == 10);
    expect_true (tx.minute () == 20);
    expect_true (tx.second () == 30);

    expect_error (Time (-1, 20, 30));
    expect_error (Time (10, -5, 30));
    expect_error (Time (10, 20, -30));
  }

  test_that ("Time as integer") {
    Time tx = Time (5000);
    expect_true (tx.seconds () == 5000);

    expect_error (Time (-1000));
  }

  test_that ("Time as string") {
    std::string t = "10:20:30";
    Time tx = Time (t);
    expect_true (tx.hour () == 10);
    expect_true (tx.minute () == 20);
    expect_true (tx.second () == 30);
    t = "1:10:20";
    expect_error (Time (t));
  }

  test_that ("UNIX timestamp to time") {
    uint64_t t = 1549764850;
    // convert it to localtime "manually" to avoid OS timezone interferring
    std::time_t tt (t);
    struct tm *t0 = localtime(&tt);

    Time tx = Time (t);
    expect_true (tx.hour () == t0->tm_hour);
    expect_true (tx.minute () == t0->tm_min);
    expect_true (tx.second () == t0->tm_sec);
  }

  test_that ("Now works") {
    Time tx = Time::now ();
    expect_true (tx.hour () >= 0);
  }

}
