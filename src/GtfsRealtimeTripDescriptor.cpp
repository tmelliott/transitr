#include <Rcpp.h>
#include <string>
#include "GtfsRealtimeTripDescriptor.h"
using namespace Rcpp;

TripDescriptor::TripDescriptor (
  Trip& trip,
  int direction_id,
  std::string start_time,
  std::string start_date,
  std::string schedule_relationship
) :
  _direction_id (direction_id),
  _start_time (start_time),
  _start_date (start_date),
  _schedule_relationship (schedule_relationship)
{
    _trip = &trip;
    _route = _trip->route();
};

Trip*
TripDescriptor::trip()
{
  return _trip;
};

Route*
TripDescriptor::route()
{
  return _route;
};

int
TripDescriptor::direction_id()
{
  return _direction_id;
};

std::string
TripDescriptor::start_time()
{
  return _start_time;
};

std::string
TripDescriptor::start_date()
{
  return _start_date;
};

std::string
TripDescriptor::schedule_relationship()
{
  return _schedule_relationship;
};

RCPP_EXPOSED_AS(Trip)

RCPP_MODULE(TripDescriptorModule) {
  class_<TripDescriptor>("TripDescriptor")

    .constructor<Trip&, int, std::string, std::string, std::string>()

    .method("trip", &TripDescriptor::trip)
    .method("route", &TripDescriptor::route)
    .method("direction_id", &TripDescriptor::direction_id)
    .method("start_time", &TripDescriptor::start_time)
    .method("start_date", &TripDescriptor::start_date)
    .method("schedule_relationship", &TripDescriptor::schedule_relationship)
    ;
}
