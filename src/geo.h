#ifndef GEO_H
#define GEO_H

#include "node.h"

#define earthRadius 6371000

struct latlng {
    double latitude;
    double longitude;
    latlng ();
    latlng (double lat, double lng);
};

double deg2rad (double deg);
double rad2deg (double rad);

double distanceEarth (double lat1d, double lon1d, double lat2d, double lon2d);
double distanceEarth (Node* from, Node* to);
double distanceEarth (latlng& p1, latlng& p2);

double bearing (double lat1d, double lon1d, double lat2d, double lon2d);
double bearing (latlng& p1, latlng& p2);

double crossTrackDistance (double latd, double lond, double latp1d, double lonp1d, double latp2d, double lonp2d);

double alongTrackDistance (double latd, double lond, double latp1d, double lonp1d, double latp2d, double lonp2d);
double alongTrackDistance (latlng& p, latlng& p1, latlng& p2);

latlng destinationPoint (latlng p0, double bearing, double distance);
Node destinationPoint (std::pair<double,double> start, double bearing, double distance);

#endif
