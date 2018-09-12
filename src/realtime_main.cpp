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

// void write_vehicles_in_parallel (Gtfs::Gtfs& gtfs, Gtfs::vehicle_map& vehicles)
// {
//     gtfs.write_vehicles (&vehicles);
//     // push sqlite -> remote postgresql
//     // {
//     //     // write to postgres in the first place (issue #5)
//     //     int rq = system ("R --slave -f scripts/copy_to_postgres.R > copy.out 2>&1 &");
//     // }
// }

// [[Rcpp::export]]
void run_realtime_model (List nw)
{
    GOOGLE_PROTOBUF_VERIFY_VERSION;

#if USE_HAVERSINE
    Rcout << "\nNOTE: using Haversine distance formula, which is slower but more accurate.\n\n";
#else
    Rcout << "\nNOTE: using Equirectangular approxmiation distance formula, which is faster but less accurate.\n\n";
#endif

    // Process nw components into c++ things
    String dbname_raw = nw["database"];
    String outputname_raw = nw["output"];
    std::string dbname (dbname_raw);
    std::string outputname (outputname_raw);
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

    // Create parameter object
    List pars = nw["parameters"];
    Gtfs::par params (pars);
    params.print ();

    // Initialize an RNG
    std::vector<RNG> rngs (params.n_core);
    unsigned int _seed = (unsigned int) time (0);
    for (int i=0; i<params.n_core; ++i) rngs.at (i).set_seed (_seed++);

    // Allow the program to be stopped gracefully    
    signal (SIGINT, intHandler);
    Timer timer;
    int tries = 0;
    while (ongoing)
    {
        timer.reset ();
        
        // call the feed once and check the result is reasonable
        int ures = rtfeed.update ();
        // 5 => "simulations completed"
        if (ures == 5) break;
        
        Rcout << "\n --- Commence iteration ---\n";
        if (ures != 0 && tries < 10)
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
        load_vehicles (&vehicles, rtfeed.feed (), &gtfs, &params);
        timer.report ("updating vehicle information");

        // Update vehicle states
        #pragma omp parallel for num_threads(params.n_core)
        for (unsigned i=0; i<vehicles.bucket_count (); ++i)
        {       
            for (auto v = vehicles.begin (i); v != vehicles.end (i); ++v)
            {
                v->second.mutate (rngs.at (omp_get_thread_num ()));
            }
        }
        timer.report ("updating vehicle states");

        // Now update the network state, using `params.n_core - 1` threads
        // std::this_thread::sleep_for (std::chrono::milliseconds (1 * 1000));
        // timer.report ("updating network state");
        
        // Predict ETAs
        #pragma omp parallel for num_threads(params.n_core)
        for (unsigned i=0; i<vehicles.bucket_count (); ++i)
        {
            for (auto v = vehicles.begin (i); v != vehicles.end (i); ++v)
            {
                v->second.predict_etas (rngs.at (omp_get_thread_num ()));
            }
        }
        timer.report ("predicting ETAs");

        // Write vehicles to (new) feed
        write_vehicles (&vehicles, outputname);
        timer.report ("writing ETAs to protobuf feed");

        gtfs.close_connection (true);
        timer.end ();

        // std::this_thread::sleep_for (std::chrono::milliseconds (10 * 1000));

    }

    Rcout << "\n\n --- Finished ---\n\n";
}
