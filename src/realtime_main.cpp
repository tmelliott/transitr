// [[Rcpp::plugins("cpp11")]]

#include <vector>
#include <string>
#include <sstream>
#include <Rcpp.h>
#include <curl/curl.h>
#include <curl/easy.h>

#include "vendor/sqlite3/sqlite3.h"

#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

#include "geo.h"
#include "gtfs.h"

using namespace Rcpp;

static size_t WriteCallback(void *contents, size_t size, size_t nmemb, void *userp)
{
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}

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
    List headers = rt["headers"];
    // std::string url (get_feed_url (nw));

    Rcout << "Oh look, it's the API url! " << url << "\n";
    CURL *curl;
    CURLcode res;
    curl_global_init (CURL_GLOBAL_DEFAULT);
    curl = curl_easy_init ();
    if (curl)
    {
        // add the headers
        struct curl_slist *chunk = NULL;

        for (int i=0; i<headers.size (); i++)
        {
            List hdr = headers[i];
            String hdrn = hdr["name"];
            String hdrv = hdr["value"];
            std::ostringstream hss;
            std::string header;
            hss << (std::string)hdrn << ": " << (std::string)hdrv;
            header = hss.str ();
            Rcout << header << "\n\n";
            chunk = curl_slist_append (chunk, header.c_str ());
        }

        res = curl_easy_setopt(curl, CURLOPT_HTTPHEADER, chunk);
        // curl_easy_setopt(curl, CURLOPT_VERBOSE, 1L);
        // 


        std::string readBuffer;
        curl_easy_setopt (curl, CURLOPT_URL, url.c_str ());
        curl_easy_setopt (curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt (curl, CURLOPT_WRITEDATA, &readBuffer);
        res = curl_easy_perform (curl);
        if (res == CURLE_OK)
        {
            Rcout << "OK!\n"
                << readBuffer << "\n";
        }
        else
        {
            Rcerr << "curl_easy_perform() failed\n";
              // curl_easy_strerror(res));
        }
        curl_easy_cleanup (curl);
        curl_slist_free_all (chunk);
    }
    curl_global_cleanup ();
    
    // Connect GTFS database
    Gtfs gtfs (dbname);




}