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
    std::vector<std::string> _urls;      // each URL will be read
    std::vector<std::string> _headers;
    transit_realtime::FeedMessage _feed;

    int _n_vehicles = 0;
    int _n_trip_updates = 0;

public:
    RealtimeFeed (std::vector<std::string>& urls, List& hdrs);

    int update ();
    transit_realtime::FeedMessage* feed ();

    int n_vehicles () { return _n_vehicles; }
    int n_trip_updates () { return _n_trip_updates; }
    
};


void load_vehicles (Gtfs::vehicle_map* vehicles,
                    transit_realtime::FeedMessage* feed,
                    Gtfs::Gtfs* gtfs, Gtfs::par* params);

void write_vehicles (Gtfs::vehicle_map* vehicles, std::string& file);


#endif
