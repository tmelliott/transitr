#include <testthat.h>
#include <string>
#include <cstdlib>

#include "realtime_feed.h"

context("GTFS classes") {
    
    // construct dummy items
    if (const char* key = std::getenv ("APIKEY"))
    {
        std::string apikey = key;
        expect_true (key != "");

        std::vector<std::string> urls;
        urls.push_back ("https://api.at.govt.nz/v2/public/realtime/vehiclelocations");
        List headers;
        {
            List keyheader = List::create (_["name"] = "Ocp-Apim-Subscription-Key",
                                           _["value"] = key);
            headers.push_back (keyheader);
        }
        {
            List typeheader = List::create (_["name"] = "Accept",
                                            _["value"] = "application/x-protobuf");
            headers.push_back (typeheader);
        }
        RealtimeFeed rt (urls, headers);

        expect_true (rt.update () == 0);
        expect_true (rt.feed ()->entity_size () > 0);
    }

}
