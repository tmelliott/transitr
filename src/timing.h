#ifndef TIMING_H
#define TIMING_H

#include <chrono>
#include <ctime>
#include <Rcpp.h>

class Timer {
private:
    std::clock_t clock;
    std::chrono::high_resolution_clock::time_point wall;

    std::clock_t clock2;
    std::chrono::high_resolution_clock::time_point wall2;

public:
    Timer ();

    void reset ();
    void report ();
    void report (const char* str);
    void end ();

    double cpu_seconds ();
    double wall_seconds ();

};

#endif
