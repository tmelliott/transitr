#ifndef GTFS_H
#define GTFS_H

#include <Rcpp.h>
#include <string>
#include <unordered_map>

#include "geo.h"
#include "time.h"

namespace Gtfs 
{

    struct ShapePt;
    struct ShapeSegment;
    struct StopTime;

    class Gtfs;
    class Agency;
    class Route;
    class Trip;
    class Shape;
    class Stop;

    struct ShapePt
    {
        latlng pt;
        double distance;
        ShapePt (double x, double y, double d);
    };
    struct ShapeSegment
    {
        // Segment* segment;
        double distance;
    };
    struct StopTime
    {
        Stop* stop;
        Trip* trip;
        Time arrival_time;
        Time departure_time;
        std::string stop_headsign;
        int pickup_type;
        int dropoff_type;
        double distance;
    };

    class Agency 
    {
    private:
        Gtfs* gtfs;
        std::string _agency_id;
        std::string _agency_name;
        std::string _agency_url;
        std::string _agency_phone;
        std::string _agency_timezone;
        std::string _agency_lang;

        bool loaded = false;    // we'll only load when requested for the first time
        bool completed = false; // set to `true` once the trip has been completed

    public:
        Agency (std::string& id, Gtfs* gtfs);

        void load ();
        void unload ();
        void unload (bool complete);

        std::string& agency_id ();
        std::string& agency_name ();
        std::string& agency_url ();
        std::string& agency_phone ();
        std::string& agency_timezone ();
        std::string& agency_lang ();
    };

    class Route 
    {
    private:
        Gtfs* gtfs;
        std::string _route_id;
        std::string _route_short_name;
        std::string _route_long_name;
        unsigned short int _route_type; // [0 tram, 1 subway/metro, 2 rail, 3 bus, 4 ferry, 5 cablecar, 6 gondola, 7 funicular]
        Agency* _agency;
        float _version;

        bool loaded = false;
        bool completed = false;

    public:
        Route (std::string& id, Gtfs* gtfs);

        void load ();
        void unload ();              // default is completed = false
        void unload (bool complete);

        std::string& route_id ();
        std::string& route_short_name ();
        std::string& route_long_name ();
        unsigned short int route_type ();
        Agency* agency ();
        float version ();
    };

    class Trip 
    {
    private:
        Gtfs* gtfs;
        std::string _trip_id;
        Route* _route;
        Shape* _shape;
        // Calendar* _calendar;
        std::vector<StopTime> _stops;
        std::string _block_id;
        bool _direction_id; // 0 or 1
        std::string _trip_headsign;
        float _version;

        bool loaded = false;
        bool completed = false;

    public:
        Trip (std::string& id, Gtfs* gtfs);

        void load ();
        void unload ();              // default is completed = false
        void unload (bool complete);

        std::string& trip_id ();
        Route* route ();
        Shape* shape ();
        // Calendar* calendar ();
        std::vector<StopTime> stops ();
        std::string& block_id ();
        bool direction_id ();
        std::string& trip_headsign ();
        float version ();
    };

    class Shape
    {
    private:
        Gtfs* gtfs;
        std::string _shape_id;
        std::vector<ShapePt> _path;
        std::vector<ShapeSegment> _segments;
        float _version;

        bool loaded = false;
        bool completed = false;

    public:
        Shape (std::string& id, Gtfs* gtfs);

        void load ();
        void unload ();
        void unload (bool complete);

        std::string& shape_id ();
        std::vector<ShapePt>& path ();
        std::vector<ShapeSegment>& segments ();
        float version ();
    };


    class Stop
    {
    private:
        Gtfs* gtfs;
        std::string _stop_id;
        latlng _stop_position;
        std::string _stop_code;
        std::string _stop_name;
        std::string _stop_desc;
        std::string _zone_id;
        std::string _parent_station;
        int _location_type;
        float _version;

        bool loaded = false;
        bool completed = false;

    public:
        Stop (std::string& id, Gtfs* gtfs);

        void load ();
        void unload ();
        void unload (bool complete);

        std::string& stop_id ();
        latlng& stop_position ();
        std::string& stop_code ();
        std::string& stop_name ();
        std::string& stop_desc ();
        std::string& zone_id ();
        std::string& parent_station ();
        int location_type ();
        float version ();
    };


    class Gtfs 
    {
    private:
        std::string _dbname;
        std::unordered_map<std::string, Agency> _agencies;
        std::unordered_map<std::string, Route> _routes;
        std::unordered_map<std::string, Trip> _trips;
        std::unordered_map<std::string, Shape> _shapes;


    public:
        Gtfs (std::string& name);

        std::string& dbname ();

        // "Find" functions
        Agency* find_agency (std::string& id);
        Route* find_route (std::string& id);
        Trip* find_trip (std::string& id);
        Shape* find_shape (std::string& id);
    };


}; // namespace Gtfs

#endif