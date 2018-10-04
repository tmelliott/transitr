#ifndef GTFS_H
#define GTFS_H

#include <string>
#include <memory>
#include <unordered_map>
#include <ctime>
#include <mutex>

#include "geo.h"
#include "time.h"
#include "rng.h"

#include "vendor/protobuf/gtfs-realtime-ext.pb.h"
#include "vendor/sqlite3/sqlite3.h"
#include <Rcpp.h>

#ifndef VERBOSE
#define VERBOSE 1
#endif

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
    class Segment;
    class Intersection;
    class Stop;
    class Calendar;

    class Vehicle;
    class Particle;

    struct eta;
    struct etaQ;
    typedef std::vector<eta> etavector;

    struct eta
    {
        std::string stop_id;
        uint64_t estimate;
        std::vector<etaQ> quantiles;
    };

    struct etaQ
    {
        float quantile;
        uint64_t time;
        etaQ (float q, uint64_t t) : quantile (q), time (t) {};
    };

    unsigned int 
    find_stop_index (double distance, std::vector<StopTime>* stops);
    unsigned int
    find_segment_index (double distance, std::vector<ShapeSegment>* segments);

    typedef std::unordered_map<std::string, Vehicle> vehicle_map;

    struct par
    {
        int n_particles = 1000;
        int n_core = 1;
        float gps_error = 5;     // std. dev. of observation error
        float system_noise = 5;  // std. dev. of speed variance/second
        float pr_stop = 0.5;
        float dwell_time = 10.0;
        float gamma = 5.0;
        bool save_timings = false;
        par () {}
        par (Rcpp::List parameters);

        void print ();
    };

    struct ShapePt
    {
        latlng pt;
        double distance;
        ShapePt (double x, double y, double d);
    };
    struct ShapeSegment
    {
        Segment* segment;
        double distance;
        ShapeSegment ();
        ShapeSegment (Segment* s, double d);
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

        StopTime ();
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

    class Segment
    {
    private:
        Gtfs* gtfs;
        int _segment_id;
        Intersection* _from;
        Intersection* _to;
        double _length;

        bool loaded = false;

        // network state
        std::vector<std::pair<int, double> > _data; // new observations as vehicles traverse network
        std::mutex data_mutex;
        uint64_t _timestamp;
        double _travel_time;
        double _uncertainty;

    public:
        Segment (int id, Gtfs* gtfs);

        void load ();
        void unload ();
        bool is_loaded ();

        int segment_id ();
        Intersection* from ();
        Intersection* to ();
        double length ();

        std::vector<std::pair<int, double> >& data ();
        uint64_t timestamp ();
        double travel_time ();
        double uncertainty ();

        double get_speed ();
        double get_speed (int delta);
        int sample_travel_time (RNG& rng);
        double sample_speed (RNG& rng);

        void push_data (int time, double err);
        std::pair<double,double> predict (int delta);
        void update (uint64_t now);
    };

    class Intersection
    {
    private:
        Gtfs* gtfs;
        int _intersection_id;
        latlng _position;

        bool loaded = false;

    public:
        Intersection (int id, Gtfs* gtfs);

        void load ();
        void unload ();

        int intersection_id ();
        latlng& position ();
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
        sqlite3* _connection = nullptr;

        std::unordered_map<std::string, Agency> _agencies;
        std::unordered_map<std::string, Route> _routes;
        std::unordered_map<std::string, Trip> _trips;
        std::unordered_map<std::string, Shape> _shapes;
        std::unordered_map<int, Segment> _segments;
        std::unordered_map<int, Intersection> _intersections;
        std::unordered_map<std::string, Stop> _stops;
        std::unordered_map<std::string, Calendar> _calendar;

    public:
        Gtfs (std::string& name);

        std::string& dbname ();
        sqlite3* get_connection ();
        void close_connection ();
        void close_connection (bool sure);
        std::unordered_map<std::string, Agency>& agencies ();
        std::unordered_map<std::string, Route>& routes ();
        std::unordered_map<std::string, Trip>& trips ();
        std::unordered_map<std::string, Shape>& shapes ();
        std::unordered_map<int, Segment>& segments ();
        std::unordered_map<int, Intersection>& intersections ();
        std::unordered_map<std::string, Stop>& stops ();
        std::unordered_map<std::string, Calendar>& calendar ();

        // "Find" functions
        Agency* find_agency (std::string& id);
        Route* find_route (std::string& id);
        Trip* find_trip (std::string& id);
        Shape* find_shape (std::string& id);
        Segment* find_segment (int id);
        Intersection* find_intersection (int id);
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
            float _gpserror;
            float _systemnoise;
            float _prstop;
            float _dwelltime;
            float _gamma;

            bool _newtrip = true;
            bool _complete = false;
            int _N;
            std::vector<Particle> _state;
            bool bad_sample;

            std::vector<unsigned int> _segment_travel_times; // segment travel times
            std::vector<uint64_t> _stop_arrival_times;     // stop arrival times
            int _current_segment;
            int _current_stop;

        public:
            Vehicle (std::string& id, par* params);

            std::string& vehicle_id ();
            Trip* trip ();
            latlng& position ();
            uint64_t timestamp ();
            unsigned delta ();

            void set_trip (Trip* trip);
            void update (const transit_realtime::VehiclePosition& vp,
                         Gtfs* gtfs);
            bool valid ();
            bool complete ();

            // statistics things
            void initialize (RNG& rng);
            void mutate (RNG& rng); // mutate state
            void select (RNG& rng); // select state (given data)
            void predict_etas (RNG& rng);
            etavector get_etas ();
            void reset ();

            double distance ();
            double speed ();
            int progress ();
            float gps_error ();
            float system_noise ();
            float pr_stop ();
            float dwell_time ();
            float gamma ();

            std::vector<unsigned int>& segment_travel_times ();
            unsigned int segment_travel_time (int l);
            int current_segment ();
            std::vector<uint64_t>& stop_arrival_times ();
            uint64_t stop_arrival_time (int m);
            int current_stop ();
    };

    class Particle {
    private:
        Vehicle* vehicle;
        double distance = 0.0;
        double speed = 0.0;
        double acceleration = 0.0;
        int accelerating = 0.0;
        std::vector<int> tt; // segment travel times
        std::vector<uint64_t> at; // stop arrival times

        bool complete = false;

        double log_likelihood = -1e6;

    public:
        Particle (double d, double s, double a, Vehicle* v);
        Particle (const Particle &p);
        ~Particle ();
        
        bool is_complete ();
        double get_distance ();
        double get_speed ();
        double get_acceleration ();
        double get_ll ();
        std::vector<uint64_t>& get_arrival_times ();
        uint64_t get_arrival_time (int i);
        std::vector<int>& get_travel_times ();
        int get_travel_time (int i);

        void travel (unsigned delta, RNG& rng);
        void predict_etas (RNG& rng);
        
        void calculate_likelihood (latlng& y, std::vector<ShapePt>* path, double sigma);
    };

}; // namespace Gtfs

#endif
