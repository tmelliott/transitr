#include "gtfs.h"
#include "vendor/sqlite3/sqlite3.h"

#include <chrono>
#include <thread>
#include <fstream>

#include "timing.h"

namespace Gtfs
{

    /***************************************************** GTFS */
    Gtfs::Gtfs (std::string& name) : _dbname (name), _startdate (time (0))
    {
        std::cout << "Connected to GTFS database `"
            << _dbname << "`\n\n"
            << " *** creating templates\n";

        // Now load all of the things ...
        sqlite3* db = get_connection ();
        if (db == nullptr)
        {
            Rcpp::Rcerr << " x Unable to connect to database\n  "
                << sqlite3_errmsg (db) << "\n";
            close_connection ();
            return;
        }

        sqlite3_stmt* stmt;
        std::string qry;
        const char* qrystr;
        // load AGENCIES
        {
            qry = "SELECT count(agency_id) FROM agency";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                close_connection ();
                return;
            }
            _agencies.reserve (sqlite3_column_int (stmt, 0));
            sqlite3_finalize (stmt);

            qry = "SELECT agency_id FROM agency";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                close_connection ();
                return;
            }
            std::string aid;
            while (sqlite3_step (stmt) == SQLITE_ROW)
            {
                aid = (char*)sqlite3_column_text (stmt, 0);
                _agencies.emplace (std::piecewise_construct,
                                   std::forward_as_tuple (aid),
                                   std::forward_as_tuple (aid, this));
            }
            sqlite3_finalize (stmt);
            std::cout << " + Created " << _agencies.size () << " agencies\n";
        }
        // load ROUTES
        {
            qry = "SELECT count(route_id) FROM routes";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                close_connection ();
                return;
            }
            _routes.reserve (sqlite3_column_int (stmt, 0));
            sqlite3_finalize (stmt);

            qry = "SELECT route_id FROM routes";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                close_connection ();
                return;
            }
            std::string rid;
            while (sqlite3_step (stmt) == SQLITE_ROW)
            {
                rid = (char*)sqlite3_column_text (stmt, 0);
                _routes.emplace (std::piecewise_construct,
                                 std::forward_as_tuple (rid),
                                 std::forward_as_tuple (rid, this));
            }
            sqlite3_finalize (stmt);
            std::cout << " + Created " << _routes.size () << " routes\n";
        }
        // load TRIPS
        {
            qry = "SELECT count(trip_id) FROM trips";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                close_connection ();
                return;
            }
            _trips.reserve (sqlite3_column_int (stmt, 0));
            sqlite3_finalize (stmt);

            // ONLY trips running today
            qry = "SELECT trip_id FROM trips";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                close_connection ();
                return;
            }
            std::string tid;
            while (sqlite3_step (stmt) == SQLITE_ROW)
            {
                tid = (char*)sqlite3_column_text (stmt, 0);
                _trips.emplace (std::piecewise_construct,
                                 std::forward_as_tuple (tid),
                                 std::forward_as_tuple (tid, this));
            }
            sqlite3_finalize (stmt);
            std::cout << " + Created " << _trips.size () << " trips\n";
        }
        // load SHAPES
        {
            qry = "SELECT count(distinct shape_id) FROM shapes";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                close_connection ();
                return;
            }
            _shapes.reserve (sqlite3_column_int (stmt, 0));
            sqlite3_finalize (stmt);

            qry ="SELECT distinct shape_id FROM shapes";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                close_connection ();
                return;
            }
            std::string sid;
            while (sqlite3_step (stmt) == SQLITE_ROW)
            {
                sid = (char*)sqlite3_column_text (stmt, 0);
                _shapes.emplace (std::piecewise_construct,
                                 std::forward_as_tuple (sid),
                                 std::forward_as_tuple (sid, this));
            }
            sqlite3_finalize (stmt);
            std::cout << " + Created " << _shapes.size () << " shapes\n";
        }
        // load SEGMENTS
        {
            qry = "SELECT count(road_segment_id) FROM road_segments";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                close_connection ();
                return;
            }
            _segments.reserve (sqlite3_column_int (stmt, 0));
            sqlite3_finalize (stmt);

            qry ="SELECT road_segment_id FROM road_segments";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                close_connection ();
                return;
            }
            int sid;
            while (sqlite3_step (stmt) == SQLITE_ROW)
            {
                sid = sqlite3_column_int (stmt, 0);
                _segments.emplace (std::piecewise_construct,
                                    std::forward_as_tuple (sid),
                                    std::forward_as_tuple (sid, this));
            }
            sqlite3_finalize (stmt);
            std::cout << " + Created " << _segments.size () << " segments\n";
        }
        // load NODES
        {
            qry = "SELECT count(node_id) FROM nodes";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                close_connection ();
                return;
            }
            _nodes.reserve (sqlite3_column_int (stmt, 0));
            sqlite3_finalize (stmt);

            qry ="SELECT node_id FROM nodes";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                close_connection ();
                return;
            }
            int iid;
            while (sqlite3_step (stmt) == SQLITE_ROW)
            {
                iid = sqlite3_column_int (stmt, 0);
                _nodes.emplace (std::piecewise_construct,
                                std::forward_as_tuple (iid),
                                std::forward_as_tuple (iid, this));
            }
            sqlite3_finalize (stmt);
            std::cout << " + Created " << _nodes.size () << " nodes\n";
        }
        // load INTERSECTIONS
        {
            qry = "SELECT count(intersection_id) FROM intersections";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                close_connection ();
                return;
            }
            _intersections.reserve (sqlite3_column_int (stmt, 0));
            sqlite3_finalize (stmt);

            qry ="SELECT intersection_id FROM intersections";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                close_connection ();
                return;
            }
            int iid;
            while (sqlite3_step (stmt) == SQLITE_ROW)
            {
                iid = sqlite3_column_int (stmt, 0);
                _intersections.emplace (std::piecewise_construct,
                                        std::forward_as_tuple (iid),
                                        std::forward_as_tuple (iid, this));
            }
            sqlite3_finalize (stmt);
            std::cout << " + Created " << _intersections.size () << " intersections\n";
        }
        // load STOPS
        {
            qry = "SELECT count(distinct stop_id) FROM stops";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                close_connection ();
                return;
            }
            _stops.reserve (sqlite3_column_int (stmt, 0));
            sqlite3_finalize (stmt);

            qry = "SELECT distinct stop_id FROM stops";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                close_connection ();
                return;
            }
            std::string sid;
            while (sqlite3_step (stmt) == SQLITE_ROW)
            {
                sid = (char*)sqlite3_column_text (stmt, 0);
                _stops.emplace (std::piecewise_construct,
                                 std::forward_as_tuple (sid),
                                 std::forward_as_tuple (sid, this));
            }
            sqlite3_finalize (stmt);
            std::cout << " + Created " << _stops.size () << " stops\n";
        }
        // load CALENDAR
        {
            qry = "SELECT count(service_id) FROM calendar";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                close_connection ();
                return;
            }
            _calendar.reserve (sqlite3_column_int (stmt, 0));
            sqlite3_finalize (stmt);

            qry = "SELECT service_id FROM calendar";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                close_connection ();
                return;
            }
            std::string cid;
            while (sqlite3_step (stmt) == SQLITE_ROW)
            {
                cid = (char*)sqlite3_column_text (stmt, 0);
                _calendar.emplace (std::piecewise_construct,
                                   std::forward_as_tuple (cid),
                                   std::forward_as_tuple (cid, this));
            }
            sqlite3_finalize (stmt);
            std::cout << " + Created " << _calendar.size () << " services\n";
        }

        close_connection ();
    }

    std::string& Gtfs::dbname ()
    {
        return _dbname;
    }

    sqlite3* Gtfs::get_connection ()
    {
        std::lock_guard<std::mutex> lk (con_lock);

        if (_connection != nullptr)
        {
            return _connection;
        }

        int maxTries = 100;
        while (maxTries > 0)
        {
            if (_connection == nullptr)
            {
                const char* dbstr = _dbname.c_str ();
                int r = sqlite3_open (dbstr, &_connection);
                if (r == SQLITE_OK)
                {
                    return _connection;
                }
            }
            maxTries--;
            std::this_thread::sleep_for (std::chrono::milliseconds (100));
        }
        std::cout << "\n max tries exceeded\n";
        return nullptr;
    }

    void Gtfs::close_connection () { close_connection (false); }
    void Gtfs::close_connection (bool sure)
    {
        if (sure)
        {
            sqlite3_close (_connection);
            _connection = nullptr;
        }
    }

    void Gtfs::set_parameters (par& params)
    {
        _parameters = params;
    }

    par* Gtfs::parameters ()
    {
        return &_parameters;
    }

    std::unordered_map<std::string, Agency>& Gtfs::agencies ()
    {
        return _agencies;
    }
    std::unordered_map<std::string, Route>& Gtfs::routes ()
    {
        return _routes;
    }
    std::unordered_map<std::string, Trip>& Gtfs::trips ()
    {
        return _trips;
    }
    std::unordered_map<std::string, Shape>& Gtfs::shapes ()
    {
        return _shapes;
    }
    std::unordered_map<int, Segment>& Gtfs::segments ()
    {
        return _segments;
    }
    std::unordered_map<int, Intersection>& Gtfs::intersections ()
    {
        return _intersections;
    }
    std::unordered_map<int, Node>& Gtfs::nodes ()
    {
        return _nodes;
    }
    std::unordered_map<std::string, Stop>& Gtfs::stops ()
    {
        return _stops;
    }
    std::unordered_map<std::string, Calendar>& Gtfs::calendar ()
    {
        return _calendar;
    }

    Agency* Gtfs::find_agency (std::string& id)
    {
        auto search = _agencies.find (id);
        if (search != _agencies.end ())
        {
            return &(search->second);
        }
        return nullptr;
    }

    Route* Gtfs::find_route (std::string& id)
    {
        auto search = _routes.find (id);
        if (search != _routes.end ())
        {
            return &(search->second);
        }
        return nullptr;
    }

    Trip* Gtfs::find_trip (std::string& id)
    {
        auto search = _trips.find (id);
        if (search != _trips.end ())
        {
            return &(search->second);
        }
        return nullptr;
    }

    Shape* Gtfs::find_shape (std::string& id)
    {
        auto search = _shapes.find (id);
        if (search != _shapes.end ())
        {
            return &(search->second);
        }
        return nullptr;
    }

    Segment* Gtfs::find_segment (int id)
    {
        auto search = _segments.find (id);
        if (search != _segments.end ())
        {
            return &(search->second);
        }
        return nullptr;
    }

    Intersection* Gtfs::find_intersection (int id)
    {
        auto search = _intersections.find (id);
        if (search != _intersections.end ())
        {
            return &(search->second);
        }
        return nullptr;
    }

    Stop* Gtfs::find_stop (std::string& id)
    {
        auto search = _stops.find (id);
        if (search != _stops.end ())
        {
            return &(search->second);
        }
        return nullptr;
    }

    Node* Gtfs::find_node (int id)
    {
        auto search = _nodes.find (id);
        if (search != _nodes.end ())
        {
            return &(search->second);
        }
        return nullptr;
    }

    Calendar* Gtfs::find_calendar (std::string& id)
    {
        auto search = _calendar.find (id);
        if (search != _calendar.end ())
        {
            return &(search->second);
        }
        return nullptr;
    }

    bool Gtfs::no_trips_remaining ()
    {
        // no easy way to do this =/
        return false;
    }


    /***************************************************** Agency */
    Agency::Agency (std::string& id, Gtfs* gtfs) :
        gtfs (gtfs), _agency_id (id)
    {
        // Rcpp::Rcout << " + Create Agency " << _agency_id << "\n";
    }

    void Agency::load ()
    {
        std::lock_guard<std::mutex> lk (load_mutex);
        if (loaded) return;
        sqlite3* db = gtfs->get_connection ();
        if (db == nullptr)
        {
            Rcpp::Rcerr << " x Unable to connect to database\n  "
                << sqlite3_errmsg (db) << "\n";
            gtfs->close_connection ();
            return;
        }

        sqlite3_stmt* stmt;
        std::string qry = "SELECT agency_name, agency_url, agency_phone, agency_timezone, agency_lang FROM agency WHERE agency_id=?";
        const char* qrystr = qry.c_str ();
        if (sqlite3_prepare_v2(db, qrystr, -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        auto astr = _agency_id.c_str ();
        if (sqlite3_bind_text (stmt, 1, astr, -1, SQLITE_STATIC) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind agency id to query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_step (stmt) != SQLITE_ROW)
        {
            Rcpp::Rcerr << " x Couldn't get agency from db\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }

        _agency_name = (char*)sqlite3_column_text (stmt, 0);
        _agency_url = (char*)sqlite3_column_text (stmt, 1);
        _agency_phone = (char*)sqlite3_column_text (stmt, 2);
        _agency_timezone = (char*)sqlite3_column_text (stmt, 3);
        _agency_lang = (char*)sqlite3_column_text (stmt, 4);

        sqlite3_finalize (stmt);
        gtfs->close_connection ();

        loaded = true;

    }

    void Agency::unload () { unload (false); }

    void Agency::unload (bool complete)
    {
        completed = complete;
        loaded = false;
        _agency_name = "";
        _agency_url = "";
        _agency_phone = "";
        _agency_timezone = "";
        _agency_lang = "";
        Rcpp::Rcout << " + Agency " << _agency_id << " is unloaded\n";
    }

    std::string& Agency::agency_id () {
        return _agency_id;
    }
    std::string& Agency::agency_name () {
        if (!loaded) load();
        return _agency_name;
    }
    std::string& Agency::agency_url () {
        if (!loaded) load();
        return _agency_url;
    }
    std::string& Agency::agency_phone () {
        if (!loaded) load();
        return _agency_phone;
    }
    std::string& Agency::agency_timezone () {
        if (!loaded) load();
        return _agency_timezone;
    }
    std::string& Agency::agency_lang () {
        if (!loaded) load();
        return _agency_lang;
    }


    /***************************************************** Route */
    Route::Route (std::string& id, Gtfs* gtfs) :
        gtfs (gtfs), _route_id (id)
    {
        // Rcpp::Rcout << " + Create Route " << _route_id << "\n";
    }

    void Route::load ()
    {
        std::lock_guard<std::mutex> lk (load_mutex);
        if (loaded) return;
        sqlite3* db = gtfs->get_connection ();
        if (db == nullptr)
        {
            Rcpp::Rcerr << " x Unable to connect to database\n  "
                << sqlite3_errmsg (db) << "\n";
            gtfs->close_connection ();
            return;
        }

        sqlite3_stmt* stmt;
        std::string qry = "SELECT route_short_name, route_long_name, route_type, agency_id, version FROM routes WHERE route_id=?";
        const char* qrystr = qry.c_str ();
        if (sqlite3_prepare_v2(db, qrystr, -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        const char* rstr = _route_id.c_str ();
        if (sqlite3_bind_text (stmt, 1, rstr, -1, SQLITE_STATIC) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind route id to query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_step (stmt) != SQLITE_ROW)
        {
            Rcpp::Rcerr << " x Couldn't get route from db\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }

        _route_short_name = (char*)sqlite3_column_text (stmt, 0);
        _route_long_name = (char*)sqlite3_column_text (stmt, 1);
        _route_type = sqlite3_column_int (stmt, 2);
        std::string agency_id = (char*)sqlite3_column_text (stmt, 3);
        _version = (float)sqlite3_column_double (stmt, 4);

        _agency = gtfs->find_agency (agency_id);

        sqlite3_finalize (stmt);
        gtfs->close_connection ();

        loaded = true;
    }

    void Route::unload () { unload (false); }

    void Route::unload (bool complete)
    {
        // only makes sense to unload things with a size
        completed = complete;
        loaded = false;
        _route_short_name = "";
        _route_long_name = "";
        _agency = nullptr;
        Rcpp::Rcout << " + Route " << _route_id << " is unloaded\n";
    }

    std::string& Route::route_id () {
        return _route_id;
    }
    std::string& Route::route_short_name () {
        if (!loaded) load();
        return _route_short_name;
    }
    std::string& Route::route_long_name () {
        if (!loaded) load();
        return _route_long_name;
    }
    unsigned short int Route::route_type () {
        if (!loaded) load();
        return _route_type;
    }
    Agency* Route::agency () {
        if (!loaded) load();
        return _agency;
    }
    float Route::version () {
        if (!loaded) load();
        return _version;
    }


    /***************************************************** Trip */
    Trip::Trip (std::string& id, Gtfs* gtfs) :
        gtfs (gtfs), _trip_id (id)
    {
    }

    void Trip::load ()
    {
        // std::lock_guard<std::recursive_mutex> lk (load_mutex);
        load_mutex.lock ();
        std::cout << " * loading trip with ID = "
            << _trip_id << "\n";
        if (loaded)
        {
            load_mutex.unlock ();
            return;
        }
        sqlite3* db = gtfs->get_connection ();
        if (db == nullptr)
        {
            Rcpp::Rcerr << " x Unable to connect to database\n  "
                << sqlite3_errmsg (db) << "\n";
            gtfs->close_connection ();
            load_mutex.unlock ();
            return;
        }
        sqlite3_stmt* stmt;
        std::string qry = "SELECT route_id, shape_id, service_id, block_id, direction_id, trip_headsign, version FROM trips WHERE trip_id=?";
        const char* qrystr = qry.c_str ();
        if (sqlite3_prepare_v2(db, qrystr, -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n"
                << "TRIP ID = " << _trip_id << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            load_mutex.unlock ();
            return;
        }
        const char* tstr = _trip_id.c_str ();
        if (sqlite3_bind_text (stmt, 1, tstr, -1, SQLITE_STATIC) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind trip id to query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            load_mutex.unlock ();
            return;
        }
        if (sqlite3_step (stmt) != SQLITE_ROW)
        {
            Rcpp::Rcerr << " x Couldn't get trip from db\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            load_mutex.unlock ();
            return;
        }

        std::string routeid = (char*)sqlite3_column_text (stmt, 0);
        _route = gtfs->find_route (routeid);
        std::string shapeid = (char*)sqlite3_column_text (stmt, 1);
        _shape = gtfs->find_shape (shapeid);
        std::string serviceid = (char*)sqlite3_column_text (stmt, 2);
        _calendar = gtfs->find_calendar (serviceid);
        if (sqlite3_column_type (stmt, 3) != SQLITE_NULL)
        {
            _block_id = (char*)sqlite3_column_text (stmt, 3);
        }
        if (sqlite3_column_type (stmt, 4) != SQLITE_NULL)
        {
            _direction_id = (bool)sqlite3_column_int (stmt, 4);
        }
        if (sqlite3_column_type (stmt, 5) != SQLITE_NULL)
        {
            _trip_headsign = (char*)sqlite3_column_text (stmt, 5);
        }
        _version = (float)sqlite3_column_double (stmt, 6);

        sqlite3_finalize (stmt);

        // Load stops
        {
            sqlite3_stmt* stmt;
            const char* qrystr;

            qry = "SELECT count(stop_id) FROM stop_times WHERE trip_id=?";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                gtfs->close_connection ();
                load_mutex.unlock ();
                return;
            }
            if (sqlite3_bind_text (stmt, 1, tstr, -1, SQLITE_STATIC) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't bind stop id to query\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                gtfs->close_connection ();
                return;
            }
            if (sqlite3_step (stmt) != SQLITE_ROW)
            {
                Rcpp::Rcerr << " x Couldn't get row count from db\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                gtfs->close_connection ();
                load_mutex.unlock ();
                return;
            }
            _stops.reserve (sqlite3_column_int (stmt, 0));
            sqlite3_finalize (stmt);

            qry = "SELECT stop_id, arrival_time, departure_time, stop_headsign, pickup_type, drop_off_type, shape_dist_traveled FROM stop_times WHERE trip_id=? ORDER BY stop_sequence";
            qrystr = qry.c_str ();
            if (sqlite3_prepare_v2(db, qrystr, -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                gtfs->close_connection ();
                load_mutex.unlock ();
                return;
            }
            if (sqlite3_bind_text (stmt, 1, tstr, -1, SQLITE_STATIC) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't bind trip id to query\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                gtfs->close_connection ();
                load_mutex.unlock ();
                return;
            }

            while (sqlite3_step (stmt) == SQLITE_ROW)
            {
                std::string stopid = (char*)sqlite3_column_text (stmt, 0);
                std::string arr = (char*)sqlite3_column_text (stmt, 1);
                std::string dep = (char*)sqlite3_column_text (stmt, 2);
                std::string head;
                if (sqlite3_column_type (stmt, 3) != SQLITE_NULL)
                {
                    head = (char*)sqlite3_column_text (stmt, 3);
                }
                _stops.emplace_back (stopid, _trip_id, arr, dep, head,
                                     sqlite3_column_int (stmt, 4),
                                     sqlite3_column_int (stmt, 5),
                                     sqlite3_column_double (stmt, 6),
                                     gtfs);
            }
            _version = (float)sqlite3_column_double (stmt, 3);

            sqlite3_finalize (stmt);
        }

        // now load stop distances ...
        std::vector<ShapeNode> nodes = _shape->nodes ();
        // for (int i=0; i<nodes.size (); i++)
        // {
        //     std::cout << "\n [" << i << "] " <<
        //         nodes.at (i).distance;
        // }
        int ni = 0;
        int si = 0;
        for (auto st = _stops.begin (); st != _stops.end (); ++st)
        {
            // find next stop node
            while (nodes.at (ni).node->node_type () == 1) ni++;
            if (nodes.at (ni).node->node_id () == st->stop->node ()->node_id ())
            {
                st->distance = nodes.at (ni).distance;
                ni++;
            }
            si++;
            // something went wrong ...
            // else if (st->stop && _shape && st->distance == 0)
            // {
            //     st->distance = _shape->distance_of (st->stop->stop_position ());
            // }
        }

        // set start time
        _start_time = _stops.at (0).departure_time;
        _arrival_times.resize (_stops.size (), 0);
        _departure_times.resize (_stops.size (), 0);

        // initialize from schedule
        _eta_state.resize (1, std::make_tuple (0, 0.0));

        _timestamp = 0;
        _stop_index = 0;
        _segment_index = 0;

        loaded = true;
        load_mutex.unlock ();
    }

    void Trip::get_db_delays ()
    {
        if (!loaded) load ();

        sqlite3* db = gtfs->get_connection ();
        if (db == nullptr)
        {
            gtfs->close_connection ();
            return;
        }

        sqlite3_stmt* stmt;
        std::string qry =
            "SELECT COUNT(type) FROM sqlite_master WHERE type='table' AND name='stop_delays'";
        const char* tblcheck = qry.c_str ();
        if (sqlite3_prepare_v2 (db, tblcheck, -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_step (stmt) != SQLITE_ROW)
        {
            Rcpp::Rcerr << " x Couldn't get segment from db\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_column_int (stmt, 0) == 0)
        {
            // segment parameters table doesn't exist
            // Rcpp::Rcout << " x No parameter table found. \n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }

        // std::cout << " - found stop delays!\n";
        sqlite3_finalize (stmt);

        // table exists, grab!
        qry = "SELECT stop_sequence, avg, sd, q50 FROM stop_delays WHERE route_id=?";
        const char* delayq = qry.c_str ();
        if (sqlite3_prepare_v2(db, delayq, -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }

        const char* tidstr = truncate_id (_trip_id).c_str ();
        if (sqlite3_bind_text (stmt, 1, tidstr, -1, SQLITE_STATIC) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Cant bind TRIP_ID to query\n   "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }

        int si;
        while (sqlite3_step (stmt) == SQLITE_ROW)
        {
            si = sqlite3_column_int (stmt, 0);
            _stops.at (si-1).average_delay = sqlite3_column_double (stmt, 1);
            _stops.at (si-1).sd_delay = sqlite3_column_double (stmt, 2);
            _stops.at (si-1).median_delay = sqlite3_column_double (stmt, 3);
        }

        sqlite3_finalize (stmt);
        gtfs->close_connection ();
    }

    void Trip::unload () { unload (false); }

    void Trip::unload (bool complete)
    {
        if (!loaded) return;
        load_mutex.lock ();

        completed = complete;
        loaded = false;
        _route = nullptr;
        _shape = nullptr;
        _calendar = nullptr;
        _stops.clear ();
        _block_id = "";
        _trip_headsign = "";
        _vehicle = nullptr;

        state_initialised = false;
#if VERBOSE > 0
        std::cout << " + Trip " << _trip_id << " is unloaded\n";
#endif
        load_mutex.unlock ();
    }

    void Trip::complete ()
    {
        completed = true;
        _vehicle = nullptr;
    }

    bool Trip::is_active (uint64_t& t)
    {
        /**
         * Returns TRUE if vehicle is not null ...
         * ACTUALLY silly because some trips start early/late
         * because bus is en-route
         */
        // if (_vehicle != nullptr) return true;

        /**
         * ... or if it is 30 mins before scheduled start
         * and less than an hour since scheduled end.
         */
        Time t0 (t);
        if (!loaded) load ();

        if (_calendar->today (t))
        {
            // t0 between (first departure - 30min, last arrival + 60min)
            if (t0 < Time (stops ().front ().departure_time.seconds () - 30*60) ||
                t0 > Time (stops ().back ().arrival_time.seconds () + 60*60))
                return false;
            // vehicle associated? it's active
            if (_vehicle != nullptr) return true;
            // otherwise only active 1-min before trip starts
            return t0 > Time (stops ().front ().departure_time.seconds () - 1*60);
        }
        return false;
    }

    bool Trip::is_active ()
    {
        if (_vehicle != nullptr) return true;
        return loaded && !completed;
    }

    std::string& Trip::trip_id () {
        return _trip_id;
    }
    Route* Trip::route () {
        if (!loaded) load ();
        return _route;
    }
    Shape* Trip::shape () {
        if (!loaded) load ();
        if (_shape == nullptr) throw std::runtime_error ("shape is null");
        return _shape;
    }
    Calendar* Trip::calendar () {
        if (!loaded) load ();
        return _calendar;
    }
    std::vector<StopTime>& Trip::stops ()
    {
        if (!loaded) load ();
        return _stops;
    }
    std::string& Trip::block_id () {
        if (!loaded) load ();
        return _block_id;
    }
    bool Trip::direction_id () {
        if (!loaded) load ();
        return _direction_id;
    }
    std::string& Trip::trip_headsign () {
        if (!loaded) load ();
        return _trip_headsign;
    }
    float Trip::version () {
        if (!loaded) load ();
        return _version;
    }

    uint64_t& Trip::timestamp ()
    {
        if (!loaded) load ();
        return _timestamp;
    }

    void Trip::set_arrival_time (int m, uint64_t t)
    {
        if (m < 0 || m >= _arrival_times.size ()) return;
        _arrival_times.at (m) = t;
    }
    void Trip::set_departure_time (int m, uint64_t t)
    {
        if (m < 0 || m >= _arrival_times.size ()) return;
        _departure_times.at (m) = t;
    }

    Time& Trip::start_time ()
    {
        if (!loaded) load ();
        return _start_time;
    }

    Vehicle* Trip::vehicle ()
    {
        return _vehicle;
    }
    void Trip::assign_vehicle (Vehicle* vehicle)
    {
        if (_vehicle != nullptr && _vehicle != vehicle)
        {
            _vehicle->remove_trip ();
        }
        _vehicle = vehicle;
        // _vehicle->set_trip (this);
    }


    /***************************************************** Stoptime */
    StopTime::StopTime ()
    {

    }
    StopTime::StopTime (std::string& stop_id, std::string& trip_id,
                        std::string& at, std::string& dt,
                        std::string& headsign,
                        int pickup, int dropoff, double dist,
                        Gtfs* gtfs)
    {
        stop = gtfs->find_stop (stop_id);
        trip = gtfs->find_trip (trip_id);
        stop->add_trip (trip);
        arrival_time = Time (at);
        departure_time = Time (dt);
        stop_headsign = headsign;
        pickup_type = pickup;
        dropoff_type = dropoff;
        distance = dist;
    }


    /***************************************************** ShapePt */
    ShapePt::ShapePt (double x, double y, double d)
    {
        pt = latlng (x, y);
        distance = d;
    }

    /***************************************************** ShapeSegment */
    ShapeSegment::ShapeSegment ()
    {

    }
    ShapeSegment::ShapeSegment (Segment* s, double d)
    {
        segment = s;
        distance = d;
    }


    /***************************************************** ShapeNode */
    ShapeNode::ShapeNode ()
    {

    }
    ShapeNode::ShapeNode (Node* n, double d)
    {
        node = n;
        distance = d;
    }


    /***************************************************** Shape */
    Shape::Shape (std::string& id, Gtfs* gtfs) :
        gtfs (gtfs), _shape_id (id) {}

    void Shape::load ()
    {
        std::lock_guard<std::mutex> lk (load_mutex);
        if (loaded) return;
        sqlite3* db = gtfs->get_connection ();
        if (db == nullptr)
        {
            Rcpp::Rcerr << " x Unable to connect to database\n  "
                << sqlite3_errmsg (db) << "\n";
            gtfs->close_connection ();
            return;
        }

        sqlite3_stmt* stmt;
        const char* qrystr;

        std::string qry = "SELECT count(shape_id) FROM shapes WHERE shape_id=?";
        qrystr = qry.c_str ();
        if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        const char* shstr = _shape_id.c_str ();
        if (sqlite3_bind_text (stmt, 1, shstr, -1, SQLITE_STATIC) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind shape id to query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_step (stmt) != SQLITE_ROW)
        {
            Rcpp::Rcerr << " x Couldn't get row count from db\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        _path.reserve (sqlite3_column_int (stmt, 0));
        sqlite3_finalize (stmt);

        qry = "SELECT shape_pt_lat, shape_pt_lon, shape_dist_traveled, version FROM shapes WHERE shape_id=? ORDER BY shape_pt_sequence";
        qrystr = qry.c_str ();
        if (sqlite3_prepare_v2(db, qrystr, -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_bind_text (stmt, 1, shstr, -1, SQLITE_STATIC) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind shape id to query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }

        double d = 0.0;
        latlng px;
        while (sqlite3_step (stmt) == SQLITE_ROW)
        {
            px.latitude = sqlite3_column_double (stmt, 0);
            px.longitude = sqlite3_column_double (stmt, 1);
            if (_path.size () > 0)
            {
                d += distanceEarth (_path.back ().pt, px);
            }
            _path.emplace_back (px.latitude, px.longitude, d);
        }
        _version = (float)sqlite3_column_double (stmt, 3);

        // and then load the shape nodes
        qry = "SELECT count(shape_id) FROM shape_nodes WHERE shape_id=?";
        qrystr = qry.c_str ();
        if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_bind_text (stmt, 1, shstr, -1, SQLITE_STATIC) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind shape id to query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_step (stmt) != SQLITE_ROW)
        {
            Rcpp::Rcerr << " x Couldn't get row count from db\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        _nodes.reserve (sqlite3_column_int (stmt, 0));
        sqlite3_finalize (stmt);

        qry = "SELECT node_id, distance_traveled FROM shape_nodes WHERE shape_id=? ORDER BY node_sequence";
        qrystr = qry.c_str ();
        if (sqlite3_prepare_v2(db, qrystr, -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_bind_text (stmt, 1, shstr, -1, SQLITE_STATIC) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind shape id to query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }

        Node* nodei;
        double di;
        while (sqlite3_step (stmt) == SQLITE_ROW)
        {
            nodei = gtfs->find_node (sqlite3_column_int (stmt, 0));
            di = sqlite3_column_double (stmt, 1);
            _nodes.emplace_back (nodei, di);
        }

        sqlite3_finalize (stmt);

        // now load the road segments
        _segments.reserve (_nodes.size () - 1);
        Segment* segi;
        qry = "SELECT road_segment_id FROM road_segments WHERE node_from=? AND node_to=?";
        qrystr = qry.c_str ();
        if (sqlite3_prepare_v2 (db, qrystr, -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        for (int i=0; i<_nodes.size ()-1; i++)
        {
            if (sqlite3_bind_int (stmt, 1, _nodes.at (i).node->node_id ()) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't bind FROM node to query\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                gtfs->close_connection ();
                return;
            }
            if (sqlite3_bind_int (stmt, 2, _nodes.at (i+1).node->node_id ()) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't bind TO node to query\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                gtfs->close_connection ();
                return;
            }
            if (sqlite3_step (stmt) != SQLITE_ROW)
            {
                Rcpp::Rcerr << " x Couldn't get row count from db\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                gtfs->close_connection ();
                return;
            }
            segi = gtfs->find_segment (sqlite3_column_int (stmt, 0));
            _segments.emplace_back (segi, _nodes.at (i).distance);
            sqlite3_reset (stmt);
        }

        _segments.reserve (sqlite3_column_int (stmt, 0));
        sqlite3_finalize (stmt);

        // sqlite3_finalize (stmt);
        gtfs->close_connection ();

        loaded = true;
    }

    void Shape::unload () { unload (false); }

    void Shape::unload (bool complete)
    {
        completed = complete;
        loaded = false;
        _path.clear ();
        _segments.clear ();
        std::cout << " + Shape " << _shape_id << " is unloaded\n";
    }

    std::string& Shape::shape_id () {
        return _shape_id;
    }

    std::vector<ShapePt>& Shape::path ()
    {
        if (!loaded) load ();
        return _path;
    }

    std::vector<ShapeSegment>& Shape::segments ()
    {
        if (!loaded) load ();
        if (_segments.size () == 0)
        {
            throw std::runtime_error ("no loaded segments");
        }
        return _segments;
    }

    std::vector<ShapeNode>& Shape::nodes ()
    {
        if (!loaded) load ();
        return _nodes;
    }

    float Shape::version () {
        if (!loaded) load ();
        return _version;
    }

    double Shape::distance_of (latlng& x)
    {
        if (!loaded) load ();
        int closest = 0;
        double dmin = 1e6;
        double di;
        for (unsigned int i=0; i<_path.size (); ++i)
        {
            di = distanceEarth (_path[i].pt, x);
            if (di < dmin)
            {
                dmin = di;
                closest = i;
            }
        }

        if (dmin < 0.0001)
        {
            return _path[closest].distance;
        }

        // closest point is the first point? must be on the first piece!
        if (closest == 0)
        {
            return alongTrackDistance(x, _path[0].pt, _path[1].pt);
        }

        // closest point is the last point? must be on the last piece!
        if (closest == _path.size () - 1)
        {
            return _path[closest - 1].distance +
                alongTrackDistance(x, _path[closest - 1].pt, _path[closest].pt);
        }

        // compute the ALONG TRACK DISTANCE for
        // segments either side of closest point
        // to determing where the point lies
        double a1, a2, d1, d2;
        a1 = alongTrackDistance(x, _path[closest - 1].pt, _path[closest].pt);
        a2 = alongTrackDistance(x, _path[closest].pt, _path[closest + 1].pt);
        d1 = _path[closest].distance - _path[closest - 1].distance;
        d2 = _path[closest + 1].distance - _path[closest].distance;

        if (a1 >= 0 && a1 < d1)
        {
            return _path[closest - 1].distance + a1;
        }
        if (a2 >= 0 && a2 < d2)
        {
            return _path[closest].distance + a2;
        }

        return 0.0;

        // pA -> pB -> pC
        // p0B is closest, but on AB or BC?
        // bool forward (true);
        // if (closest >= _path.size () - 1)
        // {
        //     forward = false;
        // }
        // else if (closest != 0)
        // {
        //     double dA = distanceEarth (_path[closest-1].pt, x);
        //     double dB = distanceEarth (_path[closest+1].pt, x);
        //     forward = dB <= dA;
        // }

        // if (forward)
        // {
        //     return _path[closest].distance +
        //         alongTrackDistance(x, _path[closest].pt, _path[closest+1].pt);
        // }
        // else
        // {
        //     return _path[closest-1].distance +
        //         alongTrackDistance(x, _path[closest-1].pt, _path[closest].pt);
        // }
    }

    latlng Shape::coordinates_of (double& d, double offset)
    {
        if (!loaded) load();

        if (d <= 0)
        {
            return _path[0].pt;
        }
        if (d >= _path.back ().distance)
        {
            return _path.back ().pt;
        }

        // creep along the shape until d < d[i]
        unsigned int i = 0;
        while (d >= _path[i+1].distance) i++;

        if (i >= _path.size ()-1)
        {
            return _path.back ().pt;
        }

        double dd = d - _path[i].distance;
        if (dd < 1)
        {
            return _path[i].pt;
        }

        // distance difference
        double b = bearing (_path[i].pt, _path[i+1].pt);

        latlng dest = destinationPoint (_path[i].pt, b, dd);

        if (offset > 0)
        {
            double b2 = fmod (b + 270., 360.);
            dest = destinationPoint (dest, b2, offset);
        }

        return dest;
    }

    latlng Shape::coordinates_of (double& d)
    {
        return coordinates_of (d, 0.0);
    }

    /***************************************************** Segment */
    Segment::Segment (int id, Gtfs* gtfs) :
        gtfs (gtfs), _segment_id (id) {}

    void Segment::load ()
    {
        std::lock_guard<std::mutex> lk (load_mutex);
        if (loaded) return;
        sqlite3* db = gtfs->get_connection ();
        if (db == nullptr)
        {
            Rcpp::Rcerr << " x Unable to connect to database\n  "
                << sqlite3_errmsg (db) << "\n";
            gtfs->close_connection ();
            return;
        }

        sqlite3_stmt* stmt;
        const char* qrystr;

        std::string qry = "SELECT node_from, node_to, length FROM road_segments WHERE road_segment_id=?";
        qrystr = qry.c_str ();
        if (sqlite3_prepare_v2(db, qrystr, -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_bind_int (stmt, 1, _segment_id) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind segment id to query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_step (stmt) != SQLITE_ROW)
        {
            Rcpp::Rcerr << " x Couldn't get segment from db\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }

        int fromid = sqlite3_column_int (stmt, 0);
        _from = gtfs->find_node (fromid);
        int toid = sqlite3_column_int (stmt, 1);
        _to = gtfs->find_node (toid);
        _length = sqlite3_column_double (stmt, 2);

        sqlite3_finalize (stmt);

        _system_noise = gtfs->parameters ()->nw_system_noise;
        _measurement_error = gtfs->parameters ()->nw_measurement_error;

        loaded = true;

        // attempt to fetch parameters from `segment_parameters` table
        qry = "SELECT COUNT(type) FROM sqlite_master WHERE type='table' AND name='segment_parameters'";
        const char* tblcheck = qry.c_str ();
        if (sqlite3_prepare_v2 (db, tblcheck, -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_step (stmt) != SQLITE_ROW)
        {
            Rcpp::Rcerr << " x Couldn't get segment from db\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_column_int (stmt, 0) == 0)
        {
            // segment parameters table doesn't exist
            // Rcpp::Rcout << " x No parameter table found. \n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }

        sqlite3_finalize (stmt);

        // it exists - go grab them values!
        qry = "SELECT q, phi, speed, speed_var, max_speed FROM segment_parameters WHERE segment_id=?";
        const char* segpars = qry.c_str ();
        if (sqlite3_prepare_v2(db, segpars, -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_bind_int (stmt, 1, _segment_id) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind segment_id to query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_step (stmt) == SQLITE_ROW)
        {
            _system_noise = sqlite3_column_double (stmt, 0);
            _state_var = pow (sqlite3_column_double (stmt, 1), 2);
            _prior_speed = sqlite3_column_double (stmt, 2);
            _prior_speed_var = fmax (1.5, sqlite3_column_double (stmt, 3));
            _max_speed = sqlite3_column_double (stmt, 4);
        }

        sqlite3_finalize (stmt);
        gtfs->close_connection ();
    }

    void Segment::unload ()
    {
        _from = nullptr;
        _to = nullptr;
        _data.clear ();
        loaded = false;
    }

    bool Segment::is_loaded ()
    {
        return loaded;
    }

    int Segment::segment_id () {
        return _segment_id;
    }

    Node* Segment::from ()
    {
        if (!loaded) load ();
        return _from;
    }
    Node* Segment::to ()
    {
        if (!loaded) load ();
        return _to;
    }
    double Segment::length ()
    {
        if (!loaded) load ();
        return _length;
    }
    double Segment::max_speed ()
    {
        if (!loaded) load ();
        return _max_speed;
    }
    double Segment::min_err ()
    {
        if (!loaded) load ();
        return _min_err;
    }
    double Segment::state_var ()
    {
        if (!loaded) load ();
        return _state_var;
    }
    double Segment::system_noise ()
    {
        if (!loaded) load ();
        return _system_noise;
    }

    /**                                          Segment state functions */
    double Segment::prior_speed ()
    {
        if (!loaded) load ();
        return _prior_speed;
    };
    double Segment::prior_speed_var ()
    {
        if (!loaded) load ();
        return _prior_speed_var;
    };

    uint64_t Segment::timestamp ()
    {
        return _timestamp;
    }
    double Segment::speed ()
    {
        if (!loaded) load ();
        return _speed;
    }
    int Segment::travel_time ()
    {
        if (!loaded) load ();
        return round (_length / _speed);
    }
    double Segment::tt_uncertainty ()
    {
        if (!loaded) load ();
        // delta method, Var(g(x)) = g'(xhat)^2 * Var(x)
        return pow (_length / _speed, 2.) * _uncertainty;
    }

    double Segment::uncertainty ()
    {
        if (!loaded) load ();
        return _uncertainty;
    }

    std::vector<std::pair<double, double> >& Segment::get_data ()
    {
        return _data;
    }

    void Segment::push_data (double speed, double err, uint64_t ts)
    {
        if (!loaded) load ();
        if (speed <= 0 || speed > _max_speed) return;
        std::lock_guard<std::mutex> lk (data_mutex);
        _data.emplace_back (speed, fmax (_min_err, err));

#if SIMULATION
        // write observation to file
        std::ostringstream fname;
        fname << "history/segment_" << _segment_id << ".csv";
        std::ofstream fout;
        fout.open (fname.str ().c_str (), std::ofstream::app);
        fout << _segment_id << "," << ts << "," << speed << "," << err << "\n";
        fout.close ();
#endif
    }


    /***************************************************** Intersection */
    Intersection::Intersection (int id, Gtfs* gtfs) :
        gtfs (gtfs), _intersection_id (id) {}

    void Intersection::load ()
    {
        std::lock_guard<std::mutex> lk (load_mutex);
        if (loaded) return;
        sqlite3* db = gtfs->get_connection ();
        if (db == nullptr)
        {
            Rcpp::Rcerr << " x Unable to connect to database\n  "
                << sqlite3_errmsg (db) << "\n";
            gtfs->close_connection ();
            return;
        }

        sqlite3_stmt* stmt;
        const char* qrystr;
        std::string qry = "SELECT intersection_lat, intersection_lon FROM intersections WHERE intersection_id=?";
        qrystr = qry.c_str ();
        if (sqlite3_prepare_v2(db, qrystr, -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_bind_int (stmt, 1, _intersection_id) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind intersection id to query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_step (stmt) != SQLITE_ROW)
        {
            Rcpp::Rcerr << " x Couldn't get intersection from db\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }

        _position = latlng (sqlite3_column_double (stmt, 0),
                            sqlite3_column_double (stmt, 1));

        sqlite3_finalize (stmt);
        gtfs->close_connection ();

        loaded = true;
    }

    void Intersection::unload ()
    {
        loaded = false;
    }

    int Intersection::intersection_id () {
        return _intersection_id;
    }

    latlng& Intersection::position ()
    {
        if (!loaded) load ();
        return _position;
    }


    /***************************************************** Stop */
    Stop::Stop (std::string& id, Gtfs* gtfs) :
        gtfs (gtfs), _stop_id (id) {}

    void Stop::load ()
    {
        std::lock_guard<std::mutex> lk (load_mutex);
        if (loaded) return;
        sqlite3* db = gtfs->get_connection ();
        if (db == nullptr)
        {
            Rcpp::Rcerr << " x Unable to connect to database\n  "
                << sqlite3_errmsg (db) << "\n";
            gtfs->close_connection ();
            return;
        }

        sqlite3_stmt* stmt;
        const char* qrystr;
        std::string qry = "SELECT stop_lat, stop_lon, stop_code, stop_name, stop_desc, zone_id, parent_station, location_type, node_id, version FROM stops WHERE stop_id=?";
        qrystr = qry.c_str ();
        if (sqlite3_prepare_v2(db, qrystr, -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        const char* ststr = _stop_id.c_str ();
        if (sqlite3_bind_text (stmt, 1, ststr, -1, SQLITE_STATIC) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind stop id to query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_step (stmt) != SQLITE_ROW)
        {
            Rcpp::Rcerr << " x Couldn't get stop from db\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }

        _stop_position = latlng (sqlite3_column_double (stmt, 0),
                                 sqlite3_column_double (stmt, 1));
        std::string serviceid = (char*)sqlite3_column_text (stmt, 2);
        if (sqlite3_column_type (stmt, 2) != SQLITE_NULL)
        {
            _stop_code = (char*)sqlite3_column_text (stmt, 2);
        }
        if (sqlite3_column_type (stmt, 3) != SQLITE_NULL)
        {
            _stop_name = (char*)sqlite3_column_text (stmt, 3);
        }
        if (sqlite3_column_type (stmt, 4) != SQLITE_NULL)
        {
            _stop_desc = (char*)sqlite3_column_text (stmt, 4);
        }
        if (sqlite3_column_type (stmt, 5) != SQLITE_NULL)
        {
            _zone_id = (char*)sqlite3_column_text (stmt, 5);
        }
        if (sqlite3_column_type (stmt, 6) != SQLITE_NULL)
        {
            _parent_station = (char*)sqlite3_column_text (stmt, 6);
        }
        if (sqlite3_column_type (stmt, 7) != SQLITE_NULL)
        {
            _location_type = sqlite3_column_int (stmt, 7);
        }
        int nodeid = sqlite3_column_int (stmt, 8);
        _node = gtfs->find_node (nodeid);
        _version = (float)sqlite3_column_double (stmt, 9);

        sqlite3_finalize (stmt);

        qry =
            "SELECT COUNT(type) FROM sqlite_master WHERE type='table' AND name='dwell_times'";
        const char* tblcheck = qry.c_str ();
        if (sqlite3_prepare_v2 (db, tblcheck, -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_step (stmt) != SQLITE_ROW)
        {
            Rcpp::Rcerr << " x Couldn't get result from db\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_column_int (stmt, 0) == 0)
        {
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }

        sqlite3_finalize (stmt);

        qry = "SELECT avg, sd, q50 FROM dwell_times WHERE stop_id=?";
        const char* delayq = qry.c_str ();
        if (sqlite3_prepare_v2(db, delayq, -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }

        const char* sidstr = truncate_id (_stop_id).c_str ();
        if (sqlite3_bind_text (stmt, 1, sidstr, -1, SQLITE_STATIC) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Cant bind STOP_ID to query\n   "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }

        if (sqlite3_step (stmt) == SQLITE_ROW)
        {
            _dwell_time = sqlite3_column_double (stmt, 0);
            _dwell_time_sd = sqlite3_column_double (stmt, 1);
            _dwell_time_median = sqlite3_column_double (stmt, 2);
        }

        sqlite3_finalize (stmt);
        gtfs->close_connection ();

        loaded = true;
    }

    void Stop::unload () { unload (false); }

    void Stop::unload (bool complete)
    {
        completed = complete;
        loaded = false;
        _stop_code = "";
        _stop_name = "";
        _stop_desc = "";
        _zone_id = "";
        _parent_station = "";
        _trips.clear ();
        std::cout << " + Stop " << _stop_id << " is unloaded\n";
    }

    std::string& Stop::stop_id () {
        return _stop_id;
    }
    latlng& Stop::stop_position ()
    {
        if (!loaded) load ();
        return _stop_position;
    }
    std::string& Stop::stop_code ()
    {
        if (!loaded) load ();
        return _stop_code;
    }
    std::string& Stop::stop_name ()
    {
        if (!loaded) load ();
        return _stop_name;
    }
    std::string& Stop::stop_desc ()
    {
        if (!loaded) load ();
        return _stop_desc;
    }
    std::string& Stop::zone_id ()
    {
        if (!loaded) load ();
        return _zone_id;
    }
    std::string& Stop::parent_station ()
    {
        if (!loaded) load ();
        return _parent_station;
    }
    int Stop::location_type ()
    {
        if (!loaded) load ();
        return _location_type;
    }
    Node* Stop::node ()
    {
        if (!loaded) load ();
        return _node;
    }
    std::vector<Trip*>& Stop::trips ()
    {
        if (!loaded) load ();
        return _trips;
    }
    float Stop::version () {
        if (!loaded) load ();
        return _version;
    }

    void Stop::add_trip (Trip* t)
    {
        // if (!loaded) load ();
        // this line leaks (??)
        // _trips.push_back (*(&t));
    }


    /***************************************************** CalendarDate */
    CalendarDate::CalendarDate (std::string& date, int type) :
    date (date), exception_type (type)
    {

    }


    /***************************************************** Calendar */
    Calendar::Calendar (std::string& id, Gtfs* gtfs) :
        gtfs (gtfs), _service_id (id)
    {
    }

    void Calendar::load ()
    {
        std::lock_guard<std::mutex> lk (load_mutex);
        if (loaded) return;
        sqlite3* db = gtfs->get_connection ();
        if (db == nullptr)
        {
            Rcpp::Rcerr << " x Unable to connect to database\n  "
                << sqlite3_errmsg (db) << "\n";
            gtfs->close_connection ();
            return;
        }

        sqlite3_stmt* stmt;
        const char* qrystr;

        std::string qry = "SELECT monday, tuesday, wednesday, thursday, friday, saturday, sunday, start_date, end_date, version FROM calendar WHERE service_id=?";
        qrystr = qry.c_str ();
        if (sqlite3_prepare_v2(db, qrystr, -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        const char* srstr = _service_id.c_str ();
        if (sqlite3_bind_text (stmt, 1, srstr, -1, SQLITE_STATIC) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind service id to query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_step (stmt) != SQLITE_ROW)
        {
            Rcpp::Rcerr << " x Couldn't get calendar from db\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }

        _monday = (bool)sqlite3_column_int (stmt, 0);
        _tuesday = (bool)sqlite3_column_int (stmt, 1);
        _wednesday = (bool)sqlite3_column_int (stmt, 2);
        _thursday = (bool)sqlite3_column_int (stmt, 3);
        _friday = (bool)sqlite3_column_int (stmt, 4);
        _saturday = (bool)sqlite3_column_int (stmt, 5);
        _sunday = (bool)sqlite3_column_int (stmt, 6);
        if (sqlite3_column_type (stmt, 7) != SQLITE_NULL)
        {
            _start_date = (char*)sqlite3_column_text (stmt, 7);
        }
        if (sqlite3_column_type (stmt, 8) != SQLITE_NULL)
        {
            _end_date = (char*)sqlite3_column_text (stmt, 8);
        }
        _version = (float)sqlite3_column_double (stmt, 9);

        sqlite3_finalize (stmt);

        // Load exceptions

        gtfs->close_connection ();

        loaded = true;
    }

    void Calendar::unload () { unload (false); }

    void Calendar::unload (bool complete)
    {
        // only makes sense to unload things with a size
        completed = complete;
        loaded = false;
        _start_date = "";
        _end_date = "";
        std::cout << " + Calendar " << _service_id << " is unloaded\n";
    }

    std::string& Calendar::service_id () {
        return _service_id;
    }
    bool Calendar::monday ()
    {
        if (!loaded) load ();
        return _monday;
    }
    bool Calendar::tuesday ()
    {
        if (!loaded) load ();
        return _tuesday;
    }
    bool Calendar::wednesday ()
    {
        if (!loaded) load ();
        return _wednesday;
    }
    bool Calendar::thursday ()
    {
        if (!loaded) load ();
        return _thursday;
    }
    bool Calendar::friday ()
    {
        if (!loaded) load ();
        return _friday;
    }
    bool Calendar::saturday ()
    {
        if (!loaded) load ();
        return _saturday;
    }
    bool Calendar::sunday ()
    {
        if (!loaded) load ();
        return _sunday;
    }
    std::string& Calendar::start_date ()
    {
        if (!loaded) load ();
        return _start_date;
    }
    std::string& Calendar::end_date ()
    {
        if (!loaded) load ();
        return _end_date;
    }
    float Calendar::version () {
        if (!loaded) load();
        return _version;
    }

    std::vector<CalendarDate*>& Calendar::exceptions ()
    {
        if (!loaded) load();
        return _exceptions;
    }

    bool Calendar::weekdays ()
    {
        if (!loaded) load();
        return _monday && _tuesday && _wednesday &&
            _thursday && _friday;
    }

    bool Calendar::today (uint64_t& t)
    {
        if (!loaded) load ();

        std::time_t tx (t);
        struct tm *t0 = localtime (&tx);
        int dow = t0->tm_wday;

        bool running =
            (_sunday & dow == 0) ||
            (_monday & dow == 1) ||
            (_tuesday && dow == 2) ||
            (_wednesday && dow == 3) ||
            (_thursday && dow == 4) ||
            (_friday && dow == 5) ||
            (_saturday && dow == 6);

        // and exceptions
        return running;
    }

    /***************************************************** Calendar */
    Node::Node (int id, Gtfs* gtfs) :
        gtfs (gtfs), _node_id (id)
    {
    }

    void Node::load ()
    {
        std::lock_guard<std::mutex> lk (load_mutex);
        if (loaded) return;
        sqlite3* db = gtfs->get_connection ();
        if (db == nullptr)
        {
            Rcpp::Rcerr << " x Unable to connect to database\n  "
                << sqlite3_errmsg (db) << "\n";
            gtfs->close_connection ();
            return;
        }

        sqlite3_stmt* stmt;
        const char* qrystr;
        std::string qry = "SELECT node_type, node_lat, node_lon FROM nodes WHERE node_id=?";
        qrystr = qry.c_str ();
        if (sqlite3_prepare_v2(db, qrystr, -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_bind_int (stmt, 1, _node_id) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind intersection id to query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_step (stmt) != SQLITE_ROW)
        {
            Rcpp::Rcerr << " x Couldn't get intersection from db\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }

        _node_type = sqlite3_column_int (stmt, 0);
        _node_position = latlng (sqlite3_column_double (stmt, 1),
                                 sqlite3_column_double (stmt, 2));

        sqlite3_finalize (stmt);
        gtfs->close_connection ();

        loaded = true;
    }

    void Node::unload (bool complete)
    {
        completed = complete;
        loaded = false;
    }

    void Node::unload () { unload (false); }

    int Node::node_id ()
    {
        return _node_id;
    }

    int Node::node_type ()
    {
        if (!loaded) load ();
        return _node_type;
    }

    latlng& Node::node_position ()
    {
        if (!loaded) load ();
        return _node_position;
    }


    /***************************************************** Parameters */
    par::par (Rcpp::List parameters)
    {
        // fetch the parameters
        Rcpp::IntegerVector nc = parameters["n_core"];
        Rcpp::IntegerVector np = parameters["n_particles"];
        Rcpp::IntegerVector noisem = parameters["noise_model"];
        Rcpp::NumericVector sigx = parameters["system_noise"];
        Rcpp::NumericVector prstop = parameters["pr_stop"];
        Rcpp::NumericVector dwell = parameters["dwell_time"];
        Rcpp::NumericVector dwell_var = parameters["dwell_time_var"];
        Rcpp::NumericVector gam = parameters["gamma"];
        Rcpp::NumericVector sigy = parameters["gps_error"];
        Rcpp::NumericVector siga = parameters["arrival_error"];
        Rcpp::NumericVector sigd = parameters["departure_error"];
        Rcpp::NumericVector signwx = parameters["nw_system_noise"];
        Rcpp::NumericVector signwy = parameters["nw_measurement_error"];
        Rcpp::IntegerVector etam = parameters["eta_model"];
        Rcpp::LogicalVector tim = parameters["save_timings"];
        Rcpp::IntegerVector resm = parameters["reset_method"];

        // set the parameters
        n_core = (int) nc[0];
        n_particles = (int) np[0];
        noise_model = (int) noisem[0];
        system_noise = (float) sigx[0];
        pr_stop = (float) prstop[0];
        dwell_time = (float) dwell[0];
        dwell_time_var = (float) dwell_var[0];
        gamma = (float) gam[0];
        gps_error = (float) sigy[0];
        arrival_error = (float) siga[0];
        departure_error = (float) sigd[0];
        nw_system_noise = (float) signwx[0];
        nw_measurement_error = (float) signwy[0];
        eta_model = (int) etam[0];
        save_timings = (bool) tim[0];
        reset_method = (int) resm[0];
    }

    void par::print ()
    {
        std::cout << "\n >>> Using the following parameters:"
            << "\n - n_core = " << n_core
            << "\n - n_particles = " << n_particles
            << "\n - system_noise = " << system_noise
            << "\n - pr_stop = " << pr_stop
            << "\n - dwell_time = " << dwell_time
            << "\n - dwell_time_var = " << dwell_time_var
            << "\n - gamma = " << gamma
            << "\n - gps_error = " << gps_error
            << "\n - arrival_error = " << arrival_error
            << "\n - departure_error = " << departure_error
            << "\n - nw_system_noise = " << nw_system_noise
            << "\n - nw_measurement_error = " << nw_measurement_error
            << "\n - save_timings = " << (save_timings ? "true" : "false")
            << "\n";
    }

    Event::Event (uint64_t ts, EventType type, std::string trip, int index) :
        timestamp (ts), type (type), trip_id (trip), position (latlng ()), stop_index (index)
    {
    }

    Event::Event (uint64_t ts, EventType type, std::string trip, int index, int delay) :
        timestamp (ts), type (type), trip_id (trip), position (latlng ()), stop_index (index), delay (delay)
    {
    }

    Event::Event (uint64_t ts, EventType type, std::string trip, latlng pos) :
        timestamp (ts), type (type), trip_id (trip), position (pos), stop_index (-1)
    {
    }

    void Event::print ()
    {
        if (type == EventType::gps)
        {
            std::cout << "position update {" <<
                position.latitude << ", " << position.longitude << "}";
        }
        else
        {
            std::cout << (type == EventType::arrival ? "arrived" : "departed")
                << " stop " << (stop_index + 1)
                << " (" << delay << "s delay)";
        }
    }
    std::string Event::type_name ()
    {
        switch(type)
        {
            case EventType::gps :
                return "gps";
            case EventType::arrival :
                return "arrival";
            case EventType::departure :
                return "departure";
        }
    }

    /***************************************************** Vehicle */
    Vehicle::Vehicle (std::string& id, par* pars) :
    _vehicle_id (id)
    {
        // set the parameters here
        _params = pars;
        _gpserror = _params->gps_error;
        _systemnoise = _params->system_noise;
        _prstop = _params->pr_stop;
        _dwelltime = _params->dwell_time;
        _dwelltimevar = _params->dwell_time_var;
        _gamma = _params->gamma;
        _arrival_error = _params->arrival_error;
        _departure_error = _params->departure_error;
        _N = _params->n_particles;
        reset_method = _params->reset_method;
#if SIMULATION
        std::ostringstream x;
        x << "history/vehicle_" << _vehicle_id << ".csv";
        store_name = x.str ();
#endif
    }

    std::string& Vehicle::vehicle_id ()
    {
        return _vehicle_id;
    }

    Trip* Vehicle::trip ()
    {
        if (! this->has_trip () ) throw std::runtime_error ("trip is null");
        return _trip;
    }
    bool Vehicle::has_trip ()
    {
        return _trip != nullptr;
    }

    latlng& Vehicle::position ()
    {
        return _position;
    }

    uint64_t Vehicle::timestamp ()
    {
        return _timestamp;
    }

    unsigned Vehicle::delta ()
    {
        return _delta;
    }

    int Vehicle::current_delay ()
    {
        return _current_delay;
    }

    void Vehicle::override_timestamp (uint64_t ts)
    {
        _timestamp = ts;
    }

    void Vehicle::add_event (Event event)
    {
        new_events.push_back (event);
    }

    std::vector<Event>& Vehicle::get_events ()
    {
         return time_events;
    }
    Event* Vehicle::latest_event ()
    {
        return _latest_event;
    }

    std::vector<STU>* Vehicle::stop_time_updates ()
    {
        return &_stop_time_updates;
    }

    void Vehicle::set_trip (Trip* trip)
    {
        if (_trip != nullptr)
        {
            // is there another vehicle assigned to that trip?
            if (_trip->vehicle () != nullptr &&
                _trip->vehicle ()->vehicle_id () != _vehicle_id)
            {
                throw 50;
            }

            // std::cout << "\n -> Transfering "
            //     << _vehicle_id << " from "
            //     << _trip->route ()->route_short_name ()
            //     << " (" << _trip->route ()->route_long_name () << ")";
            // exiting trip is "completed" (and forgets vehicle)
            _trip->complete ();
            // if (trip)
            // {
            //     std::cout << " -> " << trip->route ()->route_short_name ()
            //     << " (" << trip->route ()->route_long_name () << ")";
            // }
        }
        _trip = trip;
        _trip->assign_vehicle (this);
    }

    void Vehicle::remove_trip ()
    {
        _trip = nullptr;
    }

    void Vehicle::update (const transit_realtime::VehiclePosition& vp,
                          Gtfs* gtfs)
    {
        if (!vp.has_trip ()) return;
        if (!vp.trip ().has_trip_id ()) return;
        if (!vp.has_position ()) return;
        if (!vp.has_timestamp ()) return;

        add_event (Event (vp.timestamp (),
                          EventType::gps,
                          vp.trip ().trip_id (),
                          latlng (vp.position ().latitude (),
                                  vp.position ().longitude ())));

        return;

//         if (vp.timestamp () <= _timestamp)
//         {
//             _delta = 0;
//             return;
//         }

//         if (_trip == nullptr || _trip->trip_id () != vp.trip ().trip_id ())
//         {
//             // assign trip <--> vehicle
//             std::string tid = vp.trip ().trip_id ();
//             set_trip (gtfs->find_trip (tid));
//             _newtrip = _trip != nullptr;

//             _previous_state.clear ();
//             _previous_ts = 0;
//             _state.clear ();
//             estimated_dist = 0.0;

//             if (_trip != nullptr)
//             {
//                 _stop_time_updates.clear ();
//                 _stop_time_updates.resize (_trip->stops ().size ());
//             }
//         }

//         if (_trip == nullptr)
//         {
//             throw std::runtime_error ("Trip not found");
//         }

// #if VERBOSE == 2
//         Timer timer;
// #endif
//         _trip->route ()->load ();
// #if VERBOSE == 2
//         std::cout << " - load route (" << timer.cpu_seconds () << "ms)";
//         timer.reset ();
// #endif
//         _trip->shape ()->load ();
// #if VERBOSE == 2
//         std::cout << " - load shape (" << timer.cpu_seconds () << "ms)";
// #endif

//         _position = latlng (vp.position ().latitude (),
//                             vp.position ().longitude ());

//         double est_dist = _trip->shape ()->distance_of (_position);
//         if (!_newtrip && est_dist < estimated_dist && _previous_state.size () == _N && _previous_ts > 0)
//         {
//             // std::cout << "\n x DIE v" << _vehicle_id;
//             // If current observation est dist is LESS than previous, revert state
//             _delta = vp.timestamp () - _previous_ts;
//             _state = _previous_state;
//             _previous_state.clear ();
//             _previous_ts = 0;
//         }
//         else
//         {
//             // this state is OK - keep it
//             _delta = _timestamp == 0 ? 0 : vp.timestamp () - _timestamp;
//             _previous_state = _state;
//             _previous_ts = _timestamp;
//         }
//         estimated_dist = est_dist;
//         _timestamp = vp.timestamp ();
    }

    /**
     * Update vehicle with a GTFS TripUpdate
     * @param tu   transit_realtime::TripUpdate object
     * @param gtfs pointer to the static GTFS object
     */
    void Vehicle::update (const transit_realtime::TripUpdate& tu,
                          Gtfs* gtfs)
    {
        if (!tu.has_trip ()) return;
        if (!tu.trip ().has_trip_id ()) return;
        // this next line is outdated since Events!!
        // if (_trip != nullptr && _trip->trip_id () != tu.trip ().trip_id ()) return;
        if (!tu.has_timestamp ()) return;

        for (auto stu : tu.stop_time_update ())
        {
            if (!stu.has_stop_sequence ()) continue;
            const transit_realtime::TripUpdate::StopTimeEvent* ste;
            if (stu.has_arrival ())
            {
                ste = &(stu.arrival ());
            }
            else if (stu.has_departure ())
            {
                ste = &(stu.departure ());
            }
            else
            {
                continue;
            }
            add_event (
                Event (
                    ste->time (),
                    (stu.has_arrival () ? EventType::arrival : EventType::departure),
                    tu.trip ().trip_id (),
                    stu.stop_sequence () - 1,
                    (ste->has_delay () ? ste->delay () : 0)
                )
            );
        }


        return;

        // if (_trip == nullptr) // only if trip missing
        // {
        //     // assign trip <--> vehicle
        //     std::string tid = tu.trip ().trip_id ();
        //     set_trip (gtfs->find_trip (tid));
        //     _newtrip = _trip != nullptr;

        //     _previous_state.clear ();
        //     _previous_ts = 0;
        //     _state.clear ();
        //     estimated_dist = 0.0;

        //     if (_trip != nullptr)
        //     {
        //         _stop_time_updates.clear ();
        //         _stop_time_updates.resize (_trip->stops ().size ());
        //     }
        // }

        // make this a for loop for futureproofing
        // STU* stup;
        // auto stu = tu.stop_time_update ()[0];
        // {
        //     if (!stu.has_stop_sequence ()) return;
        //     if (_stop_time_updates.size () > 0)
        //     {
        //         _last_stop_update_index = stu.stop_sequence () - 1;
        //         stup = &(_stop_time_updates.at (_last_stop_update_index));
        //         stup->timestamp = tu.timestamp ();
        //         if (stu.has_arrival ())
        //         {
        //             auto x = stu.arrival ();
        //             if (x.has_time ())
        //                 stup->arrival_time = x.time ();
        //             if (x.has_delay ())
        //                 stup->arrival_delay = x.delay ();
        //         }
        //         if (stu.has_departure ())
        //         {
        //             auto x = stu.departure ();
        //             if (x.has_time ())
        //                 stup->departure_time = x.time ();
        //             if (x.has_delay ())
        //                 stup->departure_delay = x.delay ();
        //         }
        //     }
        // }
    }

    void Vehicle::update (Gtfs* gtfs)
    {
        if (new_events.size () == 0) return;

        // sort events, replacing any duplicated timestamps with
        // the "last" event (gps < arrival < departure)
        // - completely overwrite (trip_id too)
        std::sort (new_events.begin (), new_events.end ());
        uint64_t ts = 0;
        for (auto e : new_events)
        {
            if (time_events.size ()) ts = time_events.back ().timestamp;
            if (e.timestamp < ts) continue;
            if (e.timestamp == ts)
            {
                // update ?
                EventType prev_type = time_events.back ().type;
                if (prev_type == EventType::departure) continue;
                if (prev_type == EventType::arrival &&
                    e.type != EventType::departure) continue;

                time_events.back ().trip_id = e.trip_id;
                time_events.back ().type = e.type;
                time_events.back ().stop_index = e.stop_index;
                time_events.back ().position = e.position;
            }
            else
            {
                time_events.push_back (e);
            }
        }
        new_events.clear ();
    }

    bool Vehicle::valid ()
    {
        return _trip != nullptr && _timestamp != 0;
    }

    bool Vehicle::complete ()
    {
        return _complete;
    }

    unsigned short int Vehicle::vehicle_type ()
    {
        // somehow fetch the vehicle type from the database
        return _trip->route ()->route_type ();
    }

    bool Vehicle::is_bus ()
    {
        return vehicle_type () == 3;
    }

    bool Vehicle::is_train ()
    {
        return vehicle_type () == 2;
    }

    float Vehicle::gps_error ()
    {
        return _gpserror;
    }

    float Vehicle::system_noise ()
    {
        return _systemnoise;
    }

    float Vehicle::pr_stop ()
    {
        return _prstop;
    }
    float Vehicle::dwell_time ()
    {
        return _dwelltime;
    }
    float Vehicle::dwell_time_var ()
    {
        return _dwelltimevar;
    }
    float Vehicle::gamma ()
    {
        return _gamma;
    }
    double Vehicle::arrival_error ()
    {
        return _arrival_error;
    }
    double Vehicle::departure_error ()
    {
        return _departure_error;
    }

    std::vector<double>& Vehicle::segment_speed_avgs ()
    {
        return _segment_speed_avgs;
    }
    double Vehicle::segment_speed_avg (int l)
    {
        return _segment_speed_avgs.at (l);
    }
    int Vehicle::current_segment ()
    {
        return _current_segment;
    }
    std::vector<uint64_t>& Vehicle::stop_arrival_times ()
    {
        return _stop_arrival_times;
    }
    uint64_t Vehicle::stop_arrival_time (int m)
    {
        if (m < 0) m = 0;
        if (m >= _stop_arrival_times.size ())
            m = _stop_arrival_times.size () - 1;
        return _stop_arrival_times.at (m);
    }
    uint64_t Vehicle::stop_departure_time (int m)
    {
        if (m < 0) m = 0;
        if (m >= _stop_departure_times.size ())
            m = _stop_departure_times.size () - 1;
        return _stop_departure_times.at (m);
    }
    int Vehicle::current_stop ()
    {
        // if (_timestamp == 0) return 0;
        // int l = _current_segment;
        // std::cout << "\n l = " << l << ", len = "
        //     << _trip->shape ()->segments ().size ();
        // double d = _trip->shape ()->segments ().at (l).distance;
        // auto stops = &(_trip->stops ());
        // return find_stop_index (d, stops);
        return _current_stop;
    }

    Time& Vehicle::trip_start_time ()
    {
        return _trip->start_time ();
    }

    std::vector<Particle>* Vehicle::state ()
    {
        return &_state;
    };

    void Vehicle::revert_state ()
    {
        switch(reset_method)
        {
            case 1:
                {
                    _skip_observation = true;
                    action = "skip";
#if SIMULATION
                    this->store_state ("revert");
#endif
                    break;
                }
            case 2:
                {
                    action = "revert_state";
#if SIMULATION
                    this->store_state ("revert");
#endif
                    _state = _previous_state;
                    // set timestamp to previous obs
                    if (current_event_index > 0)
                        _timestamp = time_events.at (current_event_index - 1).timestamp;
                    else
                        _timestamp = 0;
#if SIMULATION
                    this->store_state("reverted");
#endif
                    bad_sample = false;
                    break;
                }
            case 3:
                {
                    // reset just the weights
                    break;
                }
        }

        // recalculate Neff, etc
    }

    bool Vehicle::is_errored ()
    {
        return ++error > 2;
    }

    Particle::Particle (double d, double s, double a, Vehicle* v)
    {
        vehicle = v;
        distance = d;
        speed = s;
        acceleration = a;
        stop_index = find_stop_index (d, &(v->trip ()->stops ()));
        segment_index = find_segment_index (d, &(v->trip ()->shape ()->segments ()));
        // initialize travel times to -1 and set to 0 when starting segment
        tt.resize (vehicle->trip ()->shape ()->segments ().size (), -1);
        at.resize (vehicle->trip ()->stops ().size (), 0);
        dt.resize (vehicle->trip ()->stops ().size (), 0);
        // weights initialized to ZERO
        weight = 0.0;
    }

    Particle::Particle (const Particle &p)
    {
        vehicle = p.vehicle;
        distance = p.distance;
        speed = p.speed;
        acceleration = p.acceleration;
        accelerating = p.accelerating;
        stop_index = p.stop_index;
        segment_index = p.segment_index;
        tt = p.tt;
        at = p.at;
        dt = p.dt;
        delta_ahead = p.delta_ahead;
        complete = p.complete;
        log_likelihood = p.log_likelihood;
        weight = 0.0;
    }

    Particle::~Particle ()
    {
        vehicle = nullptr;
        tt.clear ();
    }

    bool Particle::is_complete ()
    {
        return complete;
    }

    double Particle::get_distance ()
    {
        return distance;
    }

    double Particle::get_speed ()
    {
        return speed;
    }

    double Particle::get_acceleration ()
    {
        return acceleration;
    }

    unsigned int Particle::get_stop_index ()
    {
        return stop_index;
    }

    unsigned int Particle::get_segment_index ()
    {
        return segment_index;
    }

    double Particle::get_ll ()
    {
        return log_likelihood;
    }

    double Particle::get_weight ()
    {
        return weight;
    }

    std::vector<uint64_t>& Particle::get_arrival_times ()
    {
        return at;
    }

    uint64_t Particle::get_arrival_time (int i)
    {
        if (i >= at.size ()) return 0;
        return at.at (i);
    }

    void Particle::set_arrival_time (int i, uint64_t t)
    {
        if (i >= at.size ()) return;
        at.at (i) = t;
    }

    std::vector<uint64_t>& Particle::get_departure_times ()
    {
        return dt;
    }

    uint64_t Particle::get_departure_time (int i)
    {
        if (i >= dt.size ()) return 0;
        return dt.at (i);
    }

    void Particle::set_departure_time (int i, uint64_t t)
    {
        if (i >= dt.size ()) return;
        dt.at (i) = t;
    }

    std::vector<int>& Particle::get_travel_times ()
    {
        return tt;
    }
    int Particle::get_travel_time (int i)
    {
        if (i >= tt.size ()) return 0;
        return tt.at (i);
    }
    int Particle::get_travel_time_prediction (int i)
    {
        if (i >= ttpred.size ()) return 0;
        return ttpred.at (i);
    }

    /**
     * Initialize travel time in the specified segment. Currently used for initialization from
     * a trip update.
     * @param i segment index to start at
     */
    void Particle::init_travel_time (int i)
    {
        if (i >= tt.size ()) return;
        tt.at (i) = 0;
    }


    unsigned int
    find_stop_index (double distance, std::vector<StopTime>* stops)
    {
        /*
         * Return the index of the stop at which the particle LAST VISITED
         * (or is currently at).
         */
        if (distance <= 0.0) return 0;
        if (distance >= stops->back ().distance) return stops->size () - 1;
        unsigned int j = 0;
        while (stops->at (j+1).distance <= distance) j++;
        return j;
    }

    unsigned int
    find_segment_index (double distance, std::vector<ShapeSegment>* segments)
    {
        /**
         * Return the index of the segment at which the particle is now in.
         */
        if (segments->back ().distance == 0) return 0;
        if (segments->back ().distance <= distance) return segments->size () - 1;

        unsigned j;
        for (j = 0; j < segments->size () - 1; j++)
        {
            if (segments->at (j+1).distance > distance) break;
        }

        return j;
    }

    std::string truncate_id (std::string& s)
    {
        const std::string c = s;
        std::regex id_reg ("^([0-9]+)[-]");
        std::smatch match;
        if (std::regex_search (c.begin (), c.end (), match, id_reg))
        {
            // std::cout << "\nID: " << s << " -> " << match[1];
            return match[1];
        }
        return s;
    }


}; // namespace Gtfs
