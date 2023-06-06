#ifndef GTFS_ROUTE_H
#define GTFS_ROUTE_h

#include <Rcpp.h>
#include <string>

using namespace Rcpp;

// A class for a Route
class Route {
  private:
    std::string _route_id;
    std::string _route_short_name;
    std::string _route_long_name;
    int _route_type;
    std::string _agency_id;
    double _version;

  public:
    Route (
        std::string id,
        std::string short_name,
        std::string long_name,
        int type,
        std::string agency_id, // change to agency ptr in future
        double version
    );
    std::string route_id();
};

#endif
