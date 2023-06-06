#ifndef GTFS_TRIP_H
#define GTFS_TRIP_H

#include <Rcpp.h>
#include <string>
#include "GtfsScheduleRoute.h"

using namespace Rcpp;


// A class for a Trip
class Trip {
  private:
    std::string _trip_id;
    Route* _route;
    // Shape* _shape;
    // Service* _service;
    // Block* _block;
    int _direction_id;
    std::string _trip_headsign;
    double _version;

  public:
    // constructor
    Trip (
        std::string id,
        Route& route,
        int direction_id,
        std::string headsign,
        double version
    );

    // accessors
    std::string trip_id();
    Route* route();
    int direction_id();
    std::string trip_headsign();
    double version();
};


#endif
