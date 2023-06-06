#include <Rcpp.h>
#include <string>
#include "GtfsRealtimeVehicle.h"

using namespace Rcpp;

Vehicle::Vehicle (
  std::string id
) :
  vehicle_id (id)
{};

std::string
Vehicle::get_id()
{
  return vehicle_id;
};

RCPP_MODULE(VehicleModule) {
  class_<Vehicle>("Vehicle")

    .constructor<std::string>()

    .method("get_id", &Vehicle::get_id)
    ;
}
