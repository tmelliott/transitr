#ifndef GTFS_VEHICLE_H
#define GTFS_VEHICLE_H

#include <Rcpp.h>
#include <string>

using namespace Rcpp;

// A class for a vehicle
class Vehicle {
  private:
    std::string vehicle_id;

  public:
    Vehicle (std::string id);
    std::string get_id();
};

#endif
