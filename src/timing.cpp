#include "timing.h"

#include <fstream>

Timer::Timer ()
{
    reset ();
}

void Timer::reset ()
{
    clock = std::clock ();
    wall = std::chrono::high_resolution_clock::now ();   
    
    clock2 = std::clock ();
    wall2 = std::chrono::high_resolution_clock::now ();   
}

void Timer::report ()
{
    auto cend = std::clock ();
    auto wend = std::chrono::high_resolution_clock::now ();

    // elapsed time since last report
    auto ce = 1000 * (cend - clock2) / (double)(CLOCKS_PER_SEC);
    auto we = std::chrono::duration<double, std::milli> (wend - wall2).count ();

    // total time since reset
    auto ct = 1000 * (cend - clock) / (double)(CLOCKS_PER_SEC);
    auto wt = std::chrono::duration<double, std::milli> (wend - wall).count ();

    Rprintf("\n *    Time: %.2f ms (%.3f ms CPU time)",
            we, ce);
    Rprintf("\n * Elapsed: %.2f ms (%.3f ms CPU time)",
            wt, ct);
    Rcpp::Rcout << "\n *\n";

    if (save) {
        std::ofstream f;
        f.open (save_file, std::ofstream::app);
        f << (we) << "," << ce << "\n";
        f.close ();
    }

    clock2 = std::clock ();
    wall2 = std::chrono::high_resolution_clock::now (); 
}

void Timer::report (const char* str)
{
    Rcpp::Rcout << "\n *** Timings of " << str << "\n *";
    if (save) {
        std::ofstream f;
        f.open (save_file, std::ofstream::app);
        f << info << ",\"" << str << "\",";
        f.close ();
    }
    report ();
}

void Timer::end ()
{
    auto cend = std::clock ();
    auto wend = std::chrono::high_resolution_clock::now ();
    auto ct = 1000 * (cend - clock) / (double)(CLOCKS_PER_SEC);
    auto wt = std::chrono::duration<double, std::milli> (wend - wall).count ();
    Rprintf("\n *** Iteration timing: %.3f ms (%.3f ms CPU time)\n",
            wt, ct);
    reset ();
}

void Timer::save_to (const char* file, const char* infcols)
{
    save = true;
    save_file = (std::string) file;
    std::ofstream f;
    f.open (save_file);
    f << infcols << ",what,wall,cpu\n";
    f.close ();
}

void Timer::set_info (std::string inf)
{
    info = inf;
}

double Timer::cpu_seconds ()
{
    auto cend = std::clock ();
    auto ct = 1000 * (cend - clock) / (double)(CLOCKS_PER_SEC);

    return ct;
}

double Timer::wall_seconds ()
{
    auto wend = std::chrono::high_resolution_clock::now ();
    auto wt = std::chrono::duration<double, std::milli> (wend - wall).count ();

    return wt;
}
