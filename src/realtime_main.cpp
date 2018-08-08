// [[Rcpp::plugins("cpp11")]]

#include <vector>

#include "realtime_feed.h"

#include <Rcpp.h>
#include "vendor/sqlite3/sqlite3.h"

#include <chrono>
#include <thread>
typedef std::chrono::high_resolution_clock Clock;

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

#include "geo.h"
#include "gtfs.h"
#include "timing.h"
#include "rng.h"

// allow clean exit with C-c
#include <signal.h>
static volatile int ongoing = 1;
void intHandler (int dummy) {
    ongoing = 0;
}

using namespace Rcpp;

void write_vehicles_in_parallel (Gtfs::Gtfs& gtfs, Gtfs::vehicle_map& vehicles)
{
    gtfs.write_vehicles (&vehicles);
    // push sqlite -> remote postgresql
    {
        int rq = system ("R --slave -f scripts/copy_to_postgres.R > copy.out 2>&1 &");
    }
}

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
    Rcout << "\n * Running on " << numcore << " cores.";
    Rcout << "\n * Initializing " << numcore << " independent RNGs. ";
    std::vector<RNG> rngs (numcore);

    unsigned int _seed = (unsigned int) time (0);
    for (int i=0; i<numcore; ++i)
    {
        rngs.at (i).set_seed (_seed++);
        // Rcout << "\n RNG" << i << ": " 
        //     << rngs.at (i).rnorm () << ", "
        //     << rngs.at (i).rnorm () << ", "
        //     << rngs.at (i).rnorm () << ", ";
    }

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
                v->second.mutate (rngs.at (omp_get_thread_num ()));
                // printf ("\n - {%i/%i} vehicle %s: %.2f", 
                //          omp_get_thread_num ()+1, numcore,
                //          v->second.vehicle_id ().c_str (), 
                //          v->second.distance ());
            }
        }
        timer.report ("updating vehicle states");

        // Write vehicles to database on a separate thread while the network update occurs
        std::thread writev (write_vehicles_in_parallel, std::ref (gtfs), std::ref (vehicles));
        // timer.report ("inserting vehicle states into database");

        // Now update the network state
        // std::this_thread::sleep_for (std::chrono::milliseconds (1 * 1000));
        // timer.report ("updating network state");

        // Wait for vehicle writing to complete ...
        writev.join ();

        timer.end ();

        std::this_thread::sleep_for (std::chrono::milliseconds (10 * 1000));
    }
    
    Rcout << "\n\n --- Finished ---\n\n";
}
