#ifndef GTFS_H
#define GTFS_H

#include <string>
#include <memory>
#include <unordered_map>
#include <ctime>

#include "geo.h"
#include "time.h"
#include "rng.h"

#include "vendor/protobuf/gtfs-realtime.pb.h"
#include <Rcpp.h>

namespace Gtfs 
{

    struct ShapePt;
    struct ShapeSegment;
    struct StopTime;
    struct CalendarDate;

    class Gtfs;
    class Agency;
    class Route;
    class Trip;
    class Shape;
    class Stop;
    class Calendar;

    class Vehicle;
    class Particle;
    class ETA;

    unsigned int 
    find_stop_index (double distance, std::vector<StopTime>* stops);

    typedef std::unordered_map<std::string, Vehicle> vehicle_map;

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

        StopTime (std::string& stop_id, std::string& trip_id,
                  std::string& at, std::string& dt, 
                  std::string& headsign, 
                  int pickup, int dropoff, double dist,
                  Gtfs* gtfs);
    };
    struct CalendarDate
    {
        std::string date;
        int exception_type;
        CalendarDate (std::string& date, int type);
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
        Calendar* _calendar;
        std::vector<StopTime> _stops;
        std::string _block_id;
        bool _direction_id; // 0 or 1
        std::string _trip_headsign;
        float _version;

        Vehicle* _vehicle = nullptr;

        bool loaded = false;
        bool completed = false;

    public:
        Trip (std::string& id, Gtfs* gtfs);

        void load ();
        void unload ();              // default is completed = false
        void unload (bool complete);
        void complete ();

        std::string& trip_id ();
        Route* route ();
        Shape* shape ();
        Calendar* calendar ();
        std::vector<StopTime>& stops ();
        std::string& block_id ();
        bool direction_id ();
        std::string& trip_headsign ();
        float version ();

        Vehicle* vehicle ();
        void assign_vehicle (Vehicle* vehicle);
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

        double distance_of (latlng& pt);
        latlng coordinates_of (double& d);
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
        std::vector<Trip*> _trips;
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
        std::vector<Trip*>& trips ();
        float version ();

        void add_trip (Trip* t);
    };


    class Calendar 
    {
    private:
        Gtfs* gtfs;
        std::string _service_id;
        bool _monday;
        bool _tuesday;
        bool _wednesday;
        bool _thursday;
        bool _friday;
        bool _saturday;
        bool _sunday;
        std::string _start_date;
        std::string _end_date;
        float _version;
        std::vector<CalendarDate*> _exceptions;

        bool loaded = false;
        bool completed = false;

    public:
        Calendar (std::string& id, Gtfs* gtfs);

        void load ();
        void unload ();              // default is completed = false
        void unload (bool complete);

        std::string& service_id ();
        bool monday ();
        bool tuesday ();
        bool wednesday ();
        bool thursday ();
        bool friday ();
        bool saturday ();
        bool sunday ();
        std::string& start_date ();
        std::string& end_date ();
        float version ();
        std::vector<CalendarDate*>& exceptions ();

        bool weekdays ();
    };


    class Gtfs 
    {
    private:
        std::string _dbname;
        time_t _startdate;

        std::unordered_map<std::string, Agency> _agencies;
        std::unordered_map<std::string, Route> _routes;
        std::unordered_map<std::string, Trip> _trips;
        std::unordered_map<std::string, Shape> _shapes;
        std::unordered_map<std::string, Stop> _stops;
        std::unordered_map<std::string, Calendar> _calendar;

    public:
        Gtfs (std::string& name);

        std::string& dbname ();
        std::unordered_map<std::string, Agency>& agencies ();
        std::unordered_map<std::string, Route>& routes ();
        std::unordered_map<std::string, Trip>& trips ();
        std::unordered_map<std::string, Shape>& shapes ();
        std::unordered_map<std::string, Stop>& stops ();
        std::unordered_map<std::string, Calendar>& calendar ();

        // "Find" functions
        Agency* find_agency (std::string& id);
        Route* find_route (std::string& id);
        Trip* find_trip (std::string& id);
        Shape* find_shape (std::string& id);
        Stop* find_stop (std::string& id);
        Calendar* find_calendar (std::string& id);

        void write_vehicles (vehicle_map* vehicles);

        bool no_trips_remaining ();
    };

    class Vehicle {
        private:
            std::string _vehicle_id;
            Trip* _trip = nullptr;
            latlng _position;
            uint64_t _timestamp = 0;
            unsigned _delta;
            double _gpserror;

            bool _newtrip = true;
            int _N;
            std::vector<Particle> _state;

        public:
            Vehicle (std::string& id, int n, double err);

            std::string& vehicle_id ();
            Trip* trip ();
            latlng& position ();
            uint64_t timestamp ();
            unsigned delta ();

            void set_trip (Trip* trip);
            void update (const transit_realtime::VehiclePosition& vp,
                         Gtfs* gtfs);
            bool valid ();

            // statistics things
            void initialize (RNG& rng);
            void mutate (RNG& rng); // mutate state
            void select (RNG& rng); // select state (given data)
            void reset ();

            double distance ();
            double speed ();
            int progress ();
    };

    class Particle {
    private:
        Vehicle* vehicle;
        double distance = 0;
        double speed = 0;
        std::vector<unsigned int> tt;

        bool complete = false;

        double log_likelihood;

    public:
        Particle (double d, double s, Vehicle* v);
        Particle (const Particle &p);
        ~Particle ();
        
        double get_distance ();
        double get_speed ();
        double get_ll ();

        void travel (unsigned delta, RNG& rng);

        double calculate_likelihood (latlng& y, std::vector<ShapePt>* path, double sigma);
    };

}; // namespace Gtfs

#endif
