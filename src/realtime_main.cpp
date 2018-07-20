// [[Rcpp::plugins("cpp11")]]

#include <vector>
#include <Rcpp.h>

#include "vendor/sqlite3/sqlite3.h"

#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

#include "geo.h"
#include "gtfs.h"

using namespace Rcpp;

// [[Rcpp::export]]
void run_realtime_model (
    List nw, 
    int nparticles,
    int numcore)
{
    std::string dbname;
    {
        String dbname_raw = nw["database"];
        dbname = dbname_raw;
    }
    Gtfs gtfs (dbname);

    
}