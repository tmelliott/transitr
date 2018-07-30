// [[Rcpp::plugins("cpp11")]]

#include <vector>

#include "realtime_feed.h"

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
    GOOGLE_PROTOBUF_VERIFY_VERSION;

    // Process nw components into c++ things
    String dbname_raw = nw["database"];
    std::string dbname (dbname_raw);
    // std::string dbname (get_database_name (nw));
    
    // Construct the realtime feed object
    List apis = nw["apis"];
    List rt = apis["realtime"];
    String url_raw = rt["url"];
    std::string url (url_raw);
    List headers = rt["headers"];
    RealtimeFeed rtfeed (url, headers);

    // call the feed once and check the result is reasonable
    if (rtfeed.update () != 0)
    {
        throw std::invalid_argument ("Unable to fetch that url");
    }
    Rcout << "\nWe have loaded " << rtfeed.feed ()->entity_size () << " locations.\n";

    // Connect GTFS database
    Gtfs::Gtfs gtfs (dbname);

    for (auto tr : gtfs.trips ())
    {
        Gtfs::Trip* ti = &(tr.second);
        std::cout << "\n" << ti->trip_id ()
            << " starts at "
            << ti->stops ().begin ()->arrival_time;
    }

    Rcout << "\n --- Finished ---\n\n";
}
