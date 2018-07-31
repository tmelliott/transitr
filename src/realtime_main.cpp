// [[Rcpp::plugins("cpp11")]]

#include <vector>

#include "realtime_feed.h"

#include <Rcpp.h>
#include "vendor/sqlite3/sqlite3.h"

#include <chrono>
#include <thread>
typedef std::chrono::high_resolution_clock Clock;

#include "geo.h"
#include "gtfs.h"
#include "timing.h"

#include <signal.h>
static volatile int ongoing = 1;
void intHandler (int dummy) {
    ongoing = 0;
}

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

    // Connect GTFS database
    Gtfs::Gtfs gtfs (dbname);

    // Create vehicle container
    Gtfs::vehicle_map vehicles;

    // Allow the program to be stopped gracefully    
    signal (SIGINT, intHandler);
    Timer timer;
    int tries = 0;
    while (ongoing)
    {
        Rcout << "\n --- Commence iteration ---\n";
        timer.reset ();
        
        // call the feed once and check the result is reasonable
        if (rtfeed.update () != 0 && tries < 10)
        {
            Rcout << "\n x Unable to fetch URL. Trying again ...\n";
            tries++;
            std::this_thread::sleep_for (std::chrono::milliseconds (5 * 1000));
            continue;
        }
        tries = 0;
        Rcout << "\n + loaded " 
            << rtfeed.feed ()->entity_size () 
            << " vehicle positions.\n";
        timer.report ("loading vehicle positions");

        // Loading vehicle positions, assigning trips
        load_vehicles (&vehicles, rtfeed.feed (), &gtfs, nparticles);
        timer.report ("updating vehicle information");

        // Update vehicle states
        for (auto v: vehicles)
        {
            // Rcout << "\n - vehicle " << v.second.vehicle_id ();
            if (v.second.valid () && v.second.delta () > 0)
            {
                v.second.mutate ();
            }
        }
        timer.report ("loading shapes");

        timer.end ();

        std::this_thread::sleep_for (std::chrono::milliseconds (5 * 1000));
    }

    // Initialize active trips
    // for (auto tr : gtfs.trips ())
    // {
    //     Gtfs::Trip* ti = &(tr.second);
    //     std::cout << "\n" << ti->trip_id ()
    //         << " starts at "
    //         << ti->stops ().begin ()->arrival_time;
    // }

    Rcout << "\n\n --- Finished ---\n\n";
}
