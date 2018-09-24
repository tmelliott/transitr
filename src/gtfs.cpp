#include "gtfs.h"
#include "vendor/sqlite3/sqlite3.h"

#include <chrono>
#include <thread>

#include "timing.h"

namespace Gtfs 
{

    /***************************************************** GTFS */
    Gtfs::Gtfs (std::string& name) : _dbname (name), _startdate (time (0))
    {
        Rcpp::Rcout << "Connected to GTFS database `"
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
        // load AGENCIES
        {
            qry = "SELECT count(agency_id) FROM agency";
            if (sqlite3_prepare_v2 (db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
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
            if (sqlite3_prepare_v2 (db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
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
            Rcpp::Rcout << " + Created " << _agencies.size () << " agencies\n";
        }
        // load ROUTES
        {
            qry = "SELECT count(route_id) FROM routes";
            if (sqlite3_prepare_v2 (db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
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
            if (sqlite3_prepare_v2 (db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
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
            Rcpp::Rcout << " + Created " << _routes.size () << " routes\n";
        }
        // load TRIPS
        {
            qry = "SELECT count(trip_id) FROM trips";
            if (sqlite3_prepare_v2 (db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                close_connection ();
                return;
            }
            _trips.reserve (sqlite3_column_int (stmt, 0));
            sqlite3_finalize (stmt);

            qry = "SELECT trip_id FROM trips";
            if (sqlite3_prepare_v2 (db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
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
            Rcpp::Rcout << " + Created " << _trips.size () << " trips\n";
        }
        // load SHAPES
        {
            qry = "SELECT count(distinct shape_id) FROM shapes";
            if (sqlite3_prepare_v2 (db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
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
            if (sqlite3_prepare_v2 (db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
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
            Rcpp::Rcout << " + Created " << _shapes.size () << " shapes\n";
        }
        // load SEGMENTS
        {
            qry = "SELECT count(road_segment_id) FROM road_segments";
            if (sqlite3_prepare_v2 (db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
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
            if (sqlite3_prepare_v2 (db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
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
            Rcpp::Rcout << " + Created " << _segments.size () << " segments\n";
        }
        // load INTERSECTIONS
        {
            qry = "SELECT count(intersection_id) FROM intersections";
            if (sqlite3_prepare_v2 (db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
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
            if (sqlite3_prepare_v2 (db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
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
            Rcpp::Rcout << " + Created " << _intersections.size () << " intersections\n";
        }
        // load STOPS
        {
            qry = "SELECT count(distinct stop_id) FROM stops";
            if (sqlite3_prepare_v2 (db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
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
            if (sqlite3_prepare_v2 (db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
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
            Rcpp::Rcout << " + Created " << _stops.size () << " stops\n";
        }
        // load CALENDAR
        {
            qry = "SELECT count(service_id) FROM calendar";
            if (sqlite3_prepare_v2 (db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
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
            if (sqlite3_prepare_v2 (db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
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
            Rcpp::Rcout << " + Created " << _calendar.size () << " services\n";
        }

        close_connection ();    
    }

    std::string& Gtfs::dbname () 
    {
        return _dbname;
    }

    sqlite3* Gtfs::get_connection ()
    {
        if (_connection != nullptr)
        {
            return _connection;
        }
        
        int maxTries = 100;
        while (maxTries > 0)
        {
            if (_connection == nullptr)
            {
                int r = sqlite3_open (_dbname.c_str (), &_connection);
                if (r == SQLITE_OK)
                {
                    return _connection;
                }
            }
            maxTries--;
            std::this_thread::sleep_for (std::chrono::milliseconds (100));
        }
        Rcpp::Rcout << "\n max tries exceeded\n";
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
        if (sqlite3_prepare_v2(db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_bind_text (stmt, 1, _agency_id.c_str (),
                               -1, SQLITE_STATIC) != SQLITE_OK)
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
        if (sqlite3_prepare_v2(db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_bind_text (stmt, 1, _route_id.c_str (),
                               -1, SQLITE_STATIC) != SQLITE_OK)
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
        std::string qry = "SELECT route_id, shape_id, service_id, block_id, direction_id, trip_headsign, version FROM trips WHERE trip_id=?";
        if (sqlite3_prepare_v2(db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_bind_text (stmt, 1, _trip_id.c_str (),
                               -1, SQLITE_STATIC) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind trip id to query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return; 
        }
        if (sqlite3_step (stmt) != SQLITE_ROW)
        {
            Rcpp::Rcerr << " x Couldn't get trip from db\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
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
            qry = "SELECT count(stop_id) FROM stop_times WHERE trip_id=?";
            if (sqlite3_prepare_v2 (db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                gtfs->close_connection ();
                return;
            }
            if (sqlite3_bind_text (stmt, 1, _trip_id.c_str (),
                                   -1, SQLITE_STATIC) != SQLITE_OK)
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
                return;
            }
            _stops.reserve (sqlite3_column_int (stmt, 0));
            sqlite3_finalize (stmt);

            qry = "SELECT stop_id, arrival_time, departure_time, stop_headsign, pickup_type, drop_off_type, shape_dist_traveled FROM stop_times WHERE trip_id=? ORDER BY stop_sequence";
            if (sqlite3_prepare_v2(db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                gtfs->close_connection ();
                return;
            }
            if (sqlite3_bind_text (stmt, 1, _trip_id.c_str (),
                                   -1, SQLITE_STATIC) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't bind trip id to query\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                gtfs->close_connection ();
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

        gtfs->close_connection ();

        // now load stop distances ...
        for (auto st = _stops.begin (); st != _stops.end (); ++st)
        {
            if (st->stop && _shape && st->distance == 0)
            {
                st->distance = _shape->distance_of (st->stop->stop_position ());
            }
        }

        loaded = true;
    }

    void Trip::unload () { unload (false); }

    void Trip::unload (bool complete)
    {
        completed = complete;
        loaded = false;
        _route = nullptr;
        _shape = nullptr;
        _calendar = nullptr;
        _block_id = "";
        _trip_headsign = "";
        _vehicle = nullptr;
        Rcpp::Rcout << " + Trip " << _trip_id << " is unloaded\n";
    }

    void Trip::complete ()
    {
        completed = true;
        _vehicle = nullptr;
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

    Vehicle* Trip::vehicle ()
    {
        return _vehicle;
    }
    void Trip::assign_vehicle (Vehicle* vehicle)
    {
        _vehicle = vehicle;
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


    /***************************************************** Shape */
    Shape::Shape (std::string& id, Gtfs* gtfs) : 
        gtfs (gtfs), _shape_id (id) {}

    void Shape::load ()
    {
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
        std::string qry = "SELECT count(shape_id) FROM shapes WHERE shape_id=?";
        if (sqlite3_prepare_v2 (db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_bind_text (stmt, 1, _shape_id.c_str (),
                               -1, SQLITE_STATIC) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind trip id to query\n  "
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
        if (sqlite3_prepare_v2(db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_bind_text (stmt, 1, _shape_id.c_str (),
                               -1, SQLITE_STATIC) != SQLITE_OK)
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
        
        sqlite3_finalize (stmt);
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
        Rcpp::Rcout << " + Shape " << _shape_id << " is unloaded\n";
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
        return _segments;
    }

    float Shape::version () { 
        if (!loaded) load ();
        return _version; 
    }

    double Shape::distance_of (latlng& x)
    {
        if (!loaded) load ();
        int closest = 0;
        double dmin = 100000;
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
        
        if (dmin < 1)
        {
            return _path[closest].distance;
        }

        // pA -> pB -> pC
        // p0B is closest, but on AB or BC?
        bool forward (true);
        if (closest >= _path.size () - 1)
        {
            forward = false;
        }
        else if (closest != 0)
        {
            double dA = distanceEarth (_path[closest-1].pt, x);
            double dB = distanceEarth (_path[closest+1].pt, x);
            forward = dB <= dA;
        }

        if (forward)
        {
            return _path[closest].distance + 
                alongTrackDistance(x, _path[closest].pt, _path[closest+1].pt);
        }
        else
        {
            return _path[closest-1].distance +
                alongTrackDistance(x, _path[closest-1].pt, _path[closest].pt);
        }
    }

    latlng Shape::coordinates_of (double& d)
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
        if (dd < 0.01)
        {
            return _path[i].pt;
        }

        // distance difference
        double b = bearing (_path[i].pt, _path[i+1].pt);

        return destinationPoint (_path[i].pt, b, dd);       
    }


    /***************************************************** Segment */
    Segment::Segment (int id, Gtfs* gtfs) : 
        gtfs (gtfs), _segment_id (id) {}

    void Segment::load ()
    {
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
        std::string qry = "SELECT int_from, int_to, length FROM road_segments WHERE road_segment_id=?";
        if (sqlite3_prepare_v2(db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
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
        _from = gtfs->find_intersection (fromid);
        int toid = sqlite3_column_int (stmt, 1);
        _to = gtfs->find_intersection (toid);
        _length = sqlite3_column_double (stmt, 2);

        sqlite3_finalize (stmt);
        gtfs->close_connection ();

        loaded = true;
    }

    void Segment::unload ()
    {
        _from = nullptr;
        _to = nullptr;
        loaded = false;
    }

    int Segment::segment_id () {
        return _segment_id; 
    }


    /***************************************************** Intersection */
    Intersection::Intersection (int id, Gtfs* gtfs) : 
        gtfs (gtfs), _intersection_id (id) {}

    void Intersection::load ()
    {
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
        std::string qry = "SELECT intersection_lat, intersection_lon FROM intersections WHERE intersection_id=?";
        if (sqlite3_prepare_v2(db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
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
        std::string qry = "SELECT stop_lat, stop_lon, stop_code, stop_name, stop_desc, zone_id, parent_station, location_type, version FROM stops WHERE stop_id=?";
        if (sqlite3_prepare_v2(db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_bind_text (stmt, 1, _stop_id.c_str (),
                               -1, SQLITE_STATIC) != SQLITE_OK)
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
        _version = (float)sqlite3_column_double (stmt, 8);

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
        Rcpp::Rcout << " + Stop " << _stop_id << " is unloaded\n";
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
        std::string qry = "SELECT monday, tuesday, wednesday, thursday, friday, saturday, sunday, start_date, end_date, version FROM calendar WHERE service_id=?";
        if (sqlite3_prepare_v2(db, qry.c_str (), -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query `" << qry << "`\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            gtfs->close_connection ();
            return;
        }
        if (sqlite3_bind_text (stmt, 1, _service_id.c_str (),
                               -1, SQLITE_STATIC) != SQLITE_OK)
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
        Rcpp::Rcout << " + Calendar " << _service_id << " is unloaded\n";
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


    /***************************************************** Parameters */
    par::par (Rcpp::List parameters)
    {
        // fetch the parameters
        Rcpp::IntegerVector np = parameters["n_particles"];
        Rcpp::IntegerVector nc = parameters["n_core"];
        Rcpp::NumericVector sigy = parameters["gps_error"];
        Rcpp::NumericVector sigx = parameters["system_noise"];
        Rcpp::LogicalVector tim = parameters["save_timings"];

        // set the parameters
        n_particles = (int) np[0];
        n_core = (int) nc[0];
        gps_error = (float) sigy[0];
        system_noise = (float) sigx[0];
        save_timings = (bool) tim[0];
    }

    void par::print ()
    {
        std::cout << "\n >>> Using the following parameters:"
            << "\n - n_particles = " << n_particles
            << "\n - n_core = " << n_core
            << "\n - gps_error = " << gps_error
            << "\n - system_noise = " << system_noise
            << "\n - save_timings = " << (save_timings ? "true" : "false")
            << "\n";
    }

    /***************************************************** Vehicle */
    Vehicle::Vehicle (std::string& id, par* params) : 
    _vehicle_id (id)
    {
        // set the parameters here
        _gpserror = params->gps_error;
        _systemnoise = params->system_noise;
        _N = params->n_particles;
    }

    std::string& Vehicle::vehicle_id ()
    {
        return _vehicle_id;
    }

    Trip* Vehicle::trip ()
    {
        return _trip;
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

    void Vehicle::set_trip (Trip* trip)
    {
        if (_trip != nullptr)
        {
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
    }

    void Vehicle::update (const transit_realtime::VehiclePosition& vp,
                          Gtfs* gtfs)
    {
        if (!vp.has_trip ()) return;
        if (!vp.trip ().has_trip_id ()) return;
        if (!vp.has_position ()) return;
        if (!vp.has_timestamp ()) return;

        if (vp.timestamp () <= _timestamp) 
        {
            _delta = 0;
            return;
        }

        if (_trip == nullptr || _trip->trip_id () != vp.trip ().trip_id ())
        {
            // assign trip <--> vehicle
            std::string tid = vp.trip ().trip_id ();
            set_trip (gtfs->find_trip (tid));
            _newtrip = _trip != nullptr;
        }

        if (_trip == nullptr)
        {
            throw std::runtime_error ("Trip not found");
        }

#if VERBOSE == 2
        Timer timer;
#endif
        _trip->route ()->load ();
#if VERBOSE == 2
        std::cout << " - load route (" << timer.cpu_seconds () << "ms)";
        timer.reset ();
#endif
        _trip->shape ()->load ();
#if VERBOSE == 2
        std::cout << " - load shape (" << timer.cpu_seconds () << "ms)";
#endif

        _position = latlng (vp.position ().latitude (),
                            vp.position ().longitude ());

        _delta = _timestamp == 0 ? 0 : vp.timestamp () - _timestamp;
        _timestamp = vp.timestamp ();
    }

    bool Vehicle::valid ()
    {
        return (_position.latitude != 0.0 || 
                _position.longitude != 0.0) &&
            _trip != nullptr && _timestamp != 0;
    }

    bool Vehicle::complete ()
    {
        return _complete;
    }

    float Vehicle::gps_error ()
    {
        return _gpserror;
    }

    float Vehicle::system_noise ()
    {
        return _systemnoise;
    }



    Particle::Particle (double d, double s, double a, Vehicle* v)
    {
        vehicle = v;
        distance = d;
        speed = s;
        acceleration = a;
        at.resize (vehicle->trip ()->stops ().size (), 0);
    }

    Particle::Particle (const Particle &p)
    {
        vehicle = p.vehicle;
        distance = p.distance;
        speed = p.speed;
        acceleration = p.acceleration;
        accelerating = p.accelerating;
        tt = p.tt;
        at = p.at;
        complete = p.complete;
        log_likelihood = p.log_likelihood;
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

    double Particle::get_ll ()
    {
        return log_likelihood;
    }

    std::vector<uint64_t>& Particle::get_arrival_times ()
    {
        return at;
    }

    uint64_t Particle::get_arrival_time (int i)
    {
        return at.at (i);
    }


    unsigned int 
    find_stop_index (double distance, std::vector<StopTime>* stops)
    {
        /*
         * Return the index of the stop at which the particle LAST VISITED
         * (or is currently at).
         */
        if (distance <= 0) return 0;
        if (distance >= stops->back ().distance) return stops->size () - 1;
        unsigned int j = 0;
        while (stops->at (j+1).distance <= distance) j++;
        return j;
    }



}; // namespace Gtfs
