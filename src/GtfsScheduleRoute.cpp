#include <Rcpp.h>
#include <string>
#include "GtfsScheduleRoute.h"

using namespace Rcpp;


// A class for a vehicle
Route::Route (
  std::string id,
  std::string short_name,
  std::string long_name,
  int type,
  std::string agency_id, // change to agency ptr in future
  double version
) :
  _route_id (id),
  _route_short_name (short_name),
  _route_long_name (long_name),
  _route_type (type),
  _agency_id (agency_id),
  _version (version)
{};

std::string
Route::route_id()
{
  return _route_id;
};


RCPP_MODULE(RouteModule) {
  class_<Route>("Route")

    .constructor<std::string, std::string, std::string, int, std::string, double>("Route constructor")

    .method("route_id", &Route::route_id, "Get route id")
    ;
}
