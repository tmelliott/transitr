// [[Rcpp::plugins("cpp11")]]

#include <vector>
#include <Rcpp.h>
#include <curl/curl.h>

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
    // Process nw components into c++ things
    String dbname_raw = nw["database"];
    std::string dbname (dbname_raw);
    // std::string dbname (get_database_name (nw));

    List apis = nw["apis"];
    List rt = apis["realtime"];
    String url_raw = rt["url"];
    std::string url (url_raw);
    // std::string url (get_feed_url (nw));

    Rcout << "Oh look, it's the API url! " << url << "\n";
    CURL *curl;
    CURLcode res;
    curl_global_init (CURL_GLOBAL_DEFAULT);
    curl = curl_easy_init ();
    if (curl)
    {
        curl_easy_setopt (curl, CURLOPT_URL, url.c_str ());
        res = curl_easy_perform (curl);
        if (res == CURLE_OK)
        {
            Rcout << "OK!\n";
        }
        else
        {
            fprintf(stderr, "curl_easy_perform() failed: %s\n",
              curl_easy_strerror(res));
        }
        curl_easy_cleanup (curl);
    }
    curl_global_cleanup ();
    
    // Connect GTFS database
    Gtfs gtfs (dbname);




}