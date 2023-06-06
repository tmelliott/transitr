#include <Rcpp.h>
#include <string>
#include "GtfsScheduleTrip.h"

using namespace Rcpp;

Trip::Trip (
  std::string id,
  Route& route,
  int direction_id,
  std::string headsign,
  double version
) :
  _trip_id (id),
  _route (&route),
  _direction_id (direction_id),
  _trip_headsign (headsign),
  _version (version)
{
};

std::string
Trip::trip_id()
{
  return _trip_id;
};

Route*
Trip::route()
{
  return _route;
};

int
Trip::direction_id()
{
  return _direction_id;
};

std::string
Trip::trip_headsign()
{
  return _trip_headsign;
};

double
Trip::version()
{
  return _version;
};

RCPP_EXPOSED_AS(Route)

RCPP_MODULE(TripModule) {
  class_<Trip>("Trip")

    .constructor<std::string, Route&, int, std::string, double>()

    .method("trip_id", &Trip::trip_id)
    .method("route", &Trip::route)
    .method("direction_id", &Trip::direction_id)
    .method("trip_headsign", &Trip::trip_headsign)
    .method("version", &Trip::version)
    ;
}
