#ifndef REALTIME_FEED_H
#define REALTIME_FEED_H

#include <string>
#include <sstream>

#include <curl/curl.h>
#include <curl/easy.h>

#include "vendor/protobuf/gtfs-realtime.pb.h"

#include <Rcpp.h>

using namespace Rcpp;

static size_t WriteCallback(void *contents, size_t size, size_t nmemb, void *userp);

class RealtimeFeed {
private:
    std::string _url;
    std::vector<std::string> _headers;
    uint64_t _timestamp;

public:
    RealtimeFeed (std::string& url, List& hdrs);

    transit_realtime::FeedMessage get ();
    
};


#endif
