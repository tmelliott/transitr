// [[Rcpp::plugins("cpp11")]]

#include <vector>

#include "realtime_feed.h"

#include <Rcpp.h>
#include "vendor/sqlite3/sqlite3.h"

#include <chrono>
#include <thread>
#include <omp.h>
typedef std::chrono::high_resolution_clock Clock;

#include "geo.h"
#include "gtfs.h"
#include "timing.h"

// allow clean exit with C-c
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
    int numcore,
    double gpserror)
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

    // Initialize an RNG
    Rcout << "\n * Running on " << numcore << " cores.\n";

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
        load_vehicles (&vehicles, rtfeed.feed (), &gtfs, nparticles, gpserror);
        timer.report ("updating vehicle information");

        // Update vehicle states
        #pragma omp parallel for num_threads(numcore)
        for (unsigned i=0; i<vehicles.bucket_count (); ++i)
        {       
            for (auto v = vehicles.begin (i); v != vehicles.end (i); ++v)
            {
                if (!v->second.valid () || v->second.delta () == 0) continue;

                v->second.mutate ();
                Rprintf ("\n - {%i/%i} vehicle %s: %.2f", 
                         omp_get_thread_num (), numcore,
                         v->second.vehicle_id ().c_str (), 
                         v->second.distance ());
            }
        }
        timer.report ("updating vehicle states");

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
