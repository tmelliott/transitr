#ifndef GTFS_H
#define GTFS_H

#include <string>
#include <memory>
#include <unordered_map>
#include <ctime>
#include <fstream>
#include <mutex>
#include <regex>

#include "eigen3/Eigen/Dense"

#include "geo.h"
#include "time.h"
#include "rng.h"

#include "vendor/protobuf/gtfs-realtime-ext.pb.h"
#include "vendor/sqlite3/sqlite3.h"
#include <Rcpp.h>

#define RCPPTHREAD_OVERRIDE_COUT 1
#include "RcppThread.h"

#ifndef VERBOSE
#define VERBOSE 0
#endif

#ifndef SIMULATION
#define SIMULATION 0
#endif

namespace Gtfs
{

    // matrix formatting
    static const Eigen::IOFormat decimalMat (6, 0, ", ", "\n", "  [", "]");
    static const Eigen::IOFormat ColVec (6, 0, ", ", "\n", "  [", "]");
    static const Eigen::IOFormat tColVec (6, 0, " ", ", ", "", "", "  [", "]^T");
    static const Eigen::IOFormat inlineMat (6, 0, " ", ", ", "(", ")", "  [", "]");

    struct ShapePt;
    struct ShapeSegment;
    struct ShapeNode;
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
    class Node;

    class Vehicle;
    class Particle;

    struct eta;
    struct etaQ;
    typedef std::vector<eta> etavector;

    struct eta
    {
        std::string stop_id;
        uint64_t estimate = 0;
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

    std::string truncate_id (std::string& s);

    typedef std::unordered_map<std::string, Vehicle> vehicle_map;
    typedef std::unordered_map<std::string, Trip> trip_map;

    struct par
    {
        int n_core = 1;
        int n_particles = 1000;
        int noise_model = 0; // [0=once per iteration/segment/stop; 1=each second]
        float system_noise = 5;  // std. dev. of speed variance/second
        float pr_stop = 0.5;
        float dwell_time = 20.0;
        float dwell_time_var = 10.0;
        float gamma = 5.0;
        float gps_error = 5.0;   // std. dev. of observation error
        float arrival_error = 5.0;
        float departure_error = 5.0;
        float nw_system_noise = 0.001;
        float nw_measurement_error = 3.0;
        int eta_model = 0; // [0=particle filter, 1=gaussian]
        bool save_timings = false;
        int reset_method = 1;
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
    struct ShapeNode
    {
        Node* node;
        double distance;
        ShapeNode ();
        ShapeNode (Node* n, double d);
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

        /* information learnt from historical data */
        double average_delay = 0.0;
        double sd_delay = 0.0;
        double median_delay = 0.0;

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

        std::mutex load_mutex;

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

        std::mutex load_mutex;

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

        // trip state:
        Time _start_time;
        Eigen::MatrixXi _eta_matrix;
        etavector arrival_times;
        std::vector<std::tuple<uint64_t, double> > _eta_state;
        // std::vector<int> _link_times;
        // std::vector<int> _dwell_times;

        // Eigen::VectorXd B;
        // Eigen::MatrixXd E;
        // Eigen::MatrixXd Hseg;   // transform segment travel times to stop tts

        uint64_t _timestamp;
        int _delta;
        int _stop_index;
        int _segment_index;
        double _segment_progress;
        int _event_type;

        bool state_initialised = false;

        // not sure what these are for ... schedule, perhaps?
        std::vector<uint64_t> _arrival_times;
        std::vector<uint64_t> _departure_times;

        std::recursive_mutex load_mutex;

        Vehicle* _vehicle = nullptr;

        bool loaded = false;
        bool completed = false;

    public:
        Trip (std::string& id, Gtfs* gtfs);

        void load ();
        void unload ();              // default is completed = false
        void unload (bool complete);
        void complete ();
        bool is_active (uint64_t& t);
        bool is_active ();

        void get_db_delays ();

        std::string& trip_id ();
        Route* route ();
        Shape* shape ();
        Calendar* calendar ();
        std::vector<StopTime>& stops ();
        std::string& block_id ();
        bool direction_id ();
        std::string& trip_headsign ();
        float version ();
        uint64_t& timestamp ();

        void set_arrival_time (int m, uint64_t t);
        void set_departure_time (int m, uint64_t t);

        void update (uint64_t& t, RNG& rng);
        void forecast (RNG& rng);

        Time& start_time ();
        uint64_t get_eta (int i);
        // void initialize_etas (RNG& rng);
        // void update_etas (uint64_t& t, RNG& rng);
        // std::pair<Eigen::VectorXd, Eigen::MatrixXd> estimate_tt ();
        // std::pair<std::vector<Time>, std::vector<double> > calculate_etas ();
        etavector get_etas ();
        etavector& get_arrival_times ();

        void print_etas ();

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
        std::vector<ShapeNode> _nodes;
        float _version;

        std::mutex load_mutex;

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
        std::vector<ShapeNode>& nodes ();
        float version ();

        double distance_of (latlng& pt);
        latlng coordinates_of (double& d, double offset);
        latlng coordinates_of (double& d);
    };

    class Segment
    {
    private:
        Gtfs* gtfs;
        int _segment_id;
        Node* _from;
        Node* _to;
        double _length;

        /* values extracted either from the schedule or historical models */
        // double _prior_travel_time = 0.0;
        // double _prior_travel_time_var = 0.0;
        double _prior_speed = 30. / 3.6;
        double _prior_speed_var = 20.0 / pow (3.6, 2);

        std::mutex load_mutex;

        bool loaded = false;

        double _max_speed = 110.0 / 3.6; // 100kmh is the max speed of a bus (presumably)
        // double min_tt = 0.0;         // assuming vehicle traveling at max speed, this is the min time
        double _min_err = 0.5 / pow (3.6, 2); // minimum speed measurement error
        float _system_noise;            // rate of change of (mean) travel time
        float _state_var = 0.0;         /* variance between vehicles (baseline) */
        float _measurement_error;

        /**
         * type 0: just the mean travel time
         * type 1: include the first derivative of travel time (i.e., "change")
         */
        int _model_type = 0; // fixed, for now ...

        // network state
        std::vector<std::pair<double, double> > _data; // new observations as vehicles traverse network
        std::mutex data_mutex;
        uint64_t _timestamp = 0;

        double _speed;
        double _uncertainty;

    public:
        Segment (int id, Gtfs* gtfs);

        void load ();
        void unload ();
        bool is_loaded ();

        int segment_id ();
        Node* from ();
        Node* to ();
        double length ();
        double max_speed ();
        double min_err ();
        double state_var ();
        double system_noise ();

        /* State functions */
        double prior_speed ();
        double prior_speed_var ();

        uint64_t timestamp ();
        double speed ();
        double uncertainty ();

        /* Network model functions */
        std::vector<std::pair<double, double> >& get_data ();
        void push_data (double speed, double err, uint64_t ts);

        double speed (int delta);
        double sample_speed (RNG& rng);
        double sample_speed (RNG& rng, int delta);

        int travel_time ();
        double tt_uncertainty ();
        /** sample travel time from current road state
         * @param rng an RNG
         * @return an integer travel time
         */
        int sample_travel_time (RNG& rng);
        /** sample travel time for road state in delta seconds
         * @param rng an RNG
         * @param delta time, in seconds, for forecast
         * @return an integer travel time
         */
        int sample_travel_time (RNG& rng, int delta);

        std::pair<double, double> predict (int delta);
        std::pair<double, double> predict (uint64_t t);

        void update (par* params, Gtfs* gtfs);
        void update (uint64_t now);
    };

    class Intersection
    {
    private:
        Gtfs* gtfs;
        int _intersection_id;
        latlng _position;

        std::mutex load_mutex;

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
        Node* _node;
        std::vector<Trip*> _trips;
        float _version;

        std::mutex load_mutex;

        bool loaded = false;
        bool completed = false;

        double _dwell_time = 0;
        double _dwell_time_sd = 0;
        double _dwell_time_median = 0;

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
        Node* node ();
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

        std::mutex load_mutex;

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
        bool today (uint64_t& t);
    };


    class Node
    {
    private:
        Gtfs* gtfs;
        int _node_id;
        int _node_type;
        latlng _node_position;

        std::mutex load_mutex;

        bool loaded = false;
        bool completed = false;

    public:
        Node (int id, Gtfs* gtfs);

        void load ();
        void unload ();
        void unload (bool complete);

        int node_id ();
        int node_type ();
        latlng& node_position ();

    };


    class Gtfs
    {
    private:
        std::string _dbname;
        time_t _startdate;
        sqlite3* _connection = nullptr;
        par _parameters;

        std::mutex con_lock;

        std::unordered_map<std::string, Agency> _agencies;
        std::unordered_map<std::string, Route> _routes;
        std::unordered_map<std::string, Trip> _trips;
        std::unordered_map<std::string, Shape> _shapes;
        std::unordered_map<int, Segment> _segments;
        std::unordered_map<int, Intersection> _intersections;
        std::unordered_map<int, Node> _nodes;
        std::unordered_map<std::string, Stop> _stops;
        std::unordered_map<std::string, Calendar> _calendar;

    public:
        Gtfs (std::string& name);

        std::string& dbname ();
        sqlite3* get_connection ();
        void close_connection ();
        void close_connection (bool sure);
        void set_parameters (par& params);
        par* parameters ();
        std::unordered_map<std::string, Agency>& agencies ();
        std::unordered_map<std::string, Route>& routes ();
        std::unordered_map<std::string, Trip>& trips ();
        std::unordered_map<std::string, Shape>& shapes ();
        std::unordered_map<int, Segment>& segments ();
        std::unordered_map<int, Intersection>& intersections ();
        std::unordered_map<int, Node>& nodes ();
        std::unordered_map<std::string, Stop>& stops ();
        std::unordered_map<std::string, Calendar>& calendar ();

        // "Find" functions
        Agency* find_agency (std::string& id);
        Route* find_route (std::string& id);
        Trip* find_trip (std::string& id);
        Shape* find_shape (std::string& id);
        Segment* find_segment (int id);
        Intersection* find_intersection (int id);
        Node* find_node (int id);
        Stop* find_stop (std::string& id);
        Calendar* find_calendar (std::string& id);

        void write_vehicles (vehicle_map* vehicles);

        bool no_trips_remaining ();
    };

    struct STU {
        uint64_t timestamp = 0;
        uint64_t arrival_time = 0;
        int arrival_delay = 0;
        uint64_t departure_time = 0;
        int departure_delay = 0;

        bool used_arrival = false;
        bool used_departure = false;

        STU () {};
    };

    enum class EventType { gps, arrival, departure };
    struct Event {
        uint64_t timestamp;
        EventType type;
        std::string trip_id;
        latlng position;    // only for type == EventType::gps
        int stop_index;     // only for type == EventType::arrival or EventType::departure
        int delay = 0;      // for arrival/departure events
        bool used = false;  // once incorporated into likelihood, no longer use this event

        Event (uint64_t ts, EventType type, std::string trip, int index);
        Event (uint64_t ts, EventType type, std::string trip, int index, int delay);
        Event (uint64_t ts, EventType type, std::string trip, latlng pos);

        bool operator < (const Event& e) const
        {
            return (timestamp < e.timestamp);
        }

        void print ();
        std::string type_name ();
    };

    class Vehicle {
        private:
            std::mutex vehicle_lock;  // only one vehicle can do stuff at a time?
            std::string _vehicle_id;
            Trip* _trip = nullptr;
            latlng _position;
            double _velocity;
            double _bearing;
            double _odometer;
            int _stop_index;
            uint64_t _timestamp = 0;
            unsigned _delta;
            int _current_delay;  // the delay in seconds at the last reported stop

            par* _params;

            float _gpserror;
            float _systemnoise;
            float _prstop;
            float _dwelltime;
            float _dwelltimevar;
            float _gamma;

            float _arrival_error = 5.0;
            float _departure_error = 5.0;

            // reset method:
            // 1 = reset state completely
            // 2 = undo previous obs reweighting
            int reset_method = 1;

            std::vector<Event> new_events;  /** these get sorted and moved to time_events */
            std::vector<Event> time_events;
            unsigned current_event_index = 0; // almost makes `Event.used` redundant
            Event* _latest_event;

            bool _newtrip = true;
            bool _complete = false;
            int _N;
            double _Neff;
            int n_bad = 2;
            std::vector<Particle> _state;
            std::vector<Particle> _previous_state;
            uint64_t _previous_ts = 0;

            std::vector<STU> _stop_time_updates;
            int _last_stop_update_index = -1;
            bool _skip_observation = false;

            double estimated_dist = 0.0;
            double dist_to_route = 0.0;
            bool bad_sample;
            bool resample;
            int resample_count = 0;
            std::string action = "";

            std::vector<double> _segment_speed_avgs; // segment travel times
            std::vector<uint64_t> _stop_arrival_times;       // stop arrival times
            std::vector<uint64_t> _stop_departure_times;     // stop departure times
            int _current_segment;
            int _current_stop;

            // std::vector<double> _tt_state; // travel time state vector
            // std::vector<std::vector<double> > _tt_cov; // covariance matrix for travel time
            // uint64_t _tt_time; // time travel times were last updated

            int error = 0;

        public:
            Vehicle (std::string& id, par* pars);

            std::string& vehicle_id ();
            Trip* trip ();
            bool has_trip ();
            latlng& position ();
            double velocity ();
            double bearing ();
            double odometer ();
            uint64_t timestamp ();
            unsigned delta ();
            int current_delay ();

            // only for use in testing!
            void override_timestamp (uint64_t ts);

            par* params () { return _params; }

            int get_n () const { return _N; };

            void add_event (Event event);
            std::vector<Event>& get_events ();
            Event* latest_event ();

            std::vector<STU>* stop_time_updates ();

            void set_trip (Trip* trip);
            void remove_trip ();
            void update (const transit_realtime::VehiclePosition& vp,
                         Gtfs* gtfs);
            void update (const transit_realtime::TripUpdate& tu,
                         Gtfs* gtfs);
            void update (Gtfs* gtfs); // move new_events -> time_events, validate, etc
            bool valid ();
            bool complete ();

            unsigned short int vehicle_type ();
            bool is_bus ();
            bool is_train ();

            // statistics things
            void initialize (RNG& rng);
            void initialize (Event& e, RNG& rng);
            std::vector<Particle>* state ();
            void mutate (RNG& rng, Gtfs* gtfs); // mutate state
            void mutate_to (Event& e, RNG& rng); // mutate state
            void select (RNG& rng); // select state (given data)
            void predict_etas (RNG& rng);
            etavector get_etas ();
            void reset ();

            void revert_state ();

            double distance ();
            double speed ();
            int progress ();
            float gps_error ();
            float system_noise ();
            float pr_stop ();
            float dwell_time ();
            float dwell_time_var ();
            float gamma ();
            double arrival_error ();
            double departure_error ();

            std::vector<double>& segment_speed_avgs ();
            double segment_speed_avg (int l);
            int current_segment ();
            std::vector<uint64_t>& stop_arrival_times ();
            uint64_t stop_arrival_time (int m);
            uint64_t stop_departure_time (int m);
            int current_stop ();

            Time& trip_start_time (); // the time the trip started (using schedule)

            bool is_errored ();

#if SIMULATION
            std::ofstream f;
            std::string store_name = "";
            bool store_created = false;
            void store_state (std::string type);
#endif
    };

    class Particle {
    private:
        Vehicle* vehicle;
        double distance = 0.0;
        double speed = 0.0;
        double acceleration = 0.0;
        int accelerating = 0.0;
        unsigned int stop_index = 0;
        unsigned int segment_index = 0;
        std::vector<int> tt;      // segment travel times
        std::vector<int> ttpred;  // predicted travel times
        std::vector<uint64_t> at; // stop arrival times
        std::vector<uint64_t> dt; // stop departure times

        std::vector<int> eta;     // predicted etas (from now)

        int delta_ahead = 0; // seconds AHEAD of vehicle's timestamp

        bool complete = false;

        double log_likelihood = -1e6;
        double weight = 0.0;

    public:
        Particle (double d, double s, double a, Vehicle* v);
        Particle (const Particle &p);
        ~Particle ();

        bool is_complete ();
        double get_distance ();
        double get_speed ();
        double get_acceleration ();
        unsigned int get_stop_index ();
        unsigned int get_segment_index ();
        double get_ll ();
        double get_weight ();
        std::vector<uint64_t>& get_arrival_times ();
        uint64_t get_arrival_time (int i);
        void set_arrival_time (int i, uint64_t t);
        std::vector<uint64_t>& get_departure_times ();
        uint64_t get_departure_time (int i);
        std::vector<int>& get_travel_times ();
        void set_departure_time (int i, uint64_t t);
        int get_travel_time (int i);
        int get_travel_time_prediction (int i);

        void init_travel_time (int i);

        void travel (int delta, Event& e, RNG& rng);
        bool bus_stop (uint64_t time, RNG& rng);
        bool behind_event (Event& e, double delta);
        void predict_etas (RNG& rng);
        int calculate_stop_eta (double vel, int i, RNG& rng);
        // int calculate_segment_tt (double vel, int i, RNG& rng);

        void calculate_likelihood (latlng& y, std::vector<ShapePt>& path, double sigma);
        void calculate_likelihood (Event& e, double error);
        void set_weight (double w);

        // deprecated
        void calculate_arrival_likelihood (int index, uint64_t time, double error);
        void calculate_departure_likelihood (int index, uint64_t time, double error);
    };

}; // namespace Gtfs

#endif
