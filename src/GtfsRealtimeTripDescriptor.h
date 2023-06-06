#ifndef GTFS_TRIP_UPDATE_H
#define GTFS_TRIP_UPDATE_H

#include <Rcpp.h>
#include <string>
#include "GtfsScheduleTrip.h"

using namespace Rcpp;

// A class for a trip update
class TripDescriptor {
  private:
    Trip* _trip;
    Route* _route;
    int _direction_id;
    std::string _start_time;
    std::string _start_date;
    std::string _schedule_relationship;

  public:
    // constructor
    TripDescriptor (
      Trip& trip,
      int direction_id,
      std::string start_time,
      std::string start_date,
      std::string schedule_relationship
    );

    // accessors
    Trip* trip();
    Route* route();
    int direction_id();
    std::string start_time();
    std::string start_date();
    std::string schedule_relationship();
};

#endif
