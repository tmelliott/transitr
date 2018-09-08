#ifndef REALTIME_FEED_H
#define REALTIME_FEED_H

#include <string>
#include <sstream>
#include <fstream>

#include <curl/curl.h>
#include <curl/easy.h>

#include "gtfs.h"

#include <Rcpp.h>

using namespace Rcpp;

static size_t WriteCallback(void *contents, size_t size, size_t nmemb, void *userp);

class RealtimeFeed {
private:
    std::string _url;
    std::vector<std::string> _headers;
    transit_realtime::FeedMessage _feed;

public:
    RealtimeFeed (std::string& url, List& hdrs);

    int update ();
    transit_realtime::FeedMessage* feed ();
    
};


void load_vehicles (Gtfs::vehicle_map* vehicles,
                    transit_realtime::FeedMessage* feed,
                    Gtfs::Gtfs* gtfs, int n, double err);

void write_vehicles (Gtfs::vehicle_map* vehicles, std::string& file);


#endif
