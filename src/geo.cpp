// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <math.h>
#include <cmath> 
#include <vector>

#include "node.h"
#include "geo.h"

latlng::latlng ()
{
  latitude = 0.0;
  longitude = 0.0;
}
latlng::latlng (double lat, double lng)
{
  latitude = lat;
  longitude = lng;
}

double deg2rad(double deg) {
  return (deg * M_PI / 180.0);
}

double rad2deg(double rad) {
  return (rad * 180.0 / M_PI);
}

double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d) {
  double lat1r, lon1r, lat2r, lon2r, u, v, a, c;
  lat1r = deg2rad (lat1d);
  lon1r = deg2rad (lon1d);
  lat2r = deg2rad (lat2d);
  lon2r = deg2rad (lon2d);
  u = sin(deg2rad(lat2d - lat1d) / 2.0);
  v = sin(deg2rad(lon2d - lon1d) / 2.0);
  a = u * u + cos(lat1r) * cos(lat2r) * v * v;
  c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));
  return earthRadius * c;
}

double distanceEarth(Node* from, Node* to)
{
  double lat1d = std::get<1> (*from);
  double lon1d = std::get<0> (*from);
  double lat2d = std::get<1> (*to);
  double lon2d = std::get<0> (*to);
  return distanceEarth(lat1d, lon1d, lat2d, lon2d);
}

double distanceEarth (latlng& p1, latlng& p2)
{
  return distanceEarth (p1.latitude, p1.longitude, p2.latitude, p2.longitude);
}

double bearing (double lat1d, double lon1d, double lat2d, double lon2d) 
{
  double lat1r, lon1r, lat2r, lon2r, y, x, th;
  lat1r = deg2rad (lat1d);
  lon1r = deg2rad (lon1d);
  lat2r = deg2rad (lat2d);
  lon2r = deg2rad (lon2d);
  y = sin(lon2r - lon1r) * cos(lat2r);
  x = cos(lat1r) * sin(lat2r) - sin(lat1r) * cos(lat2r) * cos(lon2r - lon1r);
  th = rad2deg (atan2(y, x));
  return std::fmod (th + 360.0, 360.0);
}

double crossTrackDistance (double latd, double lond, double latp1d, double lonp1d, double latp2d, double lonp2d)
{
  double d13, t13, t12;
  d13 = distanceEarth (latp1d, lonp1d, latd, lond);
  // double d12 = distanceEarth (latp1d, lonp1d, latp2d, lonp2d);

  t13 = deg2rad (bearing (latp1d, lonp1d, latd, lond));
  t12 = deg2rad (bearing (latp1d, lonp1d, latp2d, lonp2d));
  return asin (sin (d13 / earthRadius) * sin (t13 - t12)) * earthRadius;
}

double alongTrackDistance (double latd, double lond, double latp1d, double lonp1d, double latp2d, double lonp2d)
{
  double d13, dxt;
  d13 = distanceEarth (latp1d, lonp1d, latd, lond);
  dxt = crossTrackDistance (latd, lond, latp1d, lonp1d, latp2d, lonp2d);

  return acos (cos (d13 / earthRadius) / cos (dxt / earthRadius)) * earthRadius;
}

Node destinationPoint (std::pair<double,double> start, double theta, double distance)
{
  double phi, lambda, phi2, lambda2;
  phi = deg2rad (std::get<1> (start));
  lambda = deg2rad (std::get<0> (start));

  phi2 = asin (sin (phi) * cos (distance / earthRadius) +
               cos (phi) * sin (distance / earthRadius) * cos (deg2rad (theta)));
  lambda2 = lambda + atan2 (sin (deg2rad (theta)) * sin (distance / earthRadius) * cos (phi),
                            cos (distance / earthRadius) - sin (phi) * sin (phi2));

  lambda2 = fmod(rad2deg (lambda2) + 540, 360) - 180;
  return Node (lambda2, rad2deg (phi2));
}
