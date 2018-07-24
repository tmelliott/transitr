#include "gtfs.h"
#include "vendor/sqlite3/sqlite3.h"

namespace Gtfs 
{

    /***************************************************** GTFS */
    Gtfs::Gtfs (std::string& name) : _dbname (name) 
    {
        Rcpp::Rcout << "Connected to GTFS database `"
            << _dbname << "`\n\n"
            << " *** creating templates\n";

        // Now load all of the things ...
        sqlite3* db;
        sqlite3_stmt* stmt;
        if (sqlite3_open (_dbname.c_str (), &db))
        {
            Rcpp::Rcerr << " x Unable to connect to database\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_close (db);
            return;
        }
        // load AGENCIES
        {
            if (sqlite3_prepare_v2 (db, "SELECT count(agency_id) FROM agency",
                                    -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                sqlite3_close (db);
                return;
            }
            _agencies.reserve (sqlite3_column_int (stmt, 0));
            sqlite3_finalize (stmt);

            if (sqlite3_prepare_v2 (db, "SELECT agency_id FROM agency",
                                    -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                sqlite3_close (db);
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
            if (sqlite3_prepare_v2 (db, "SELECT count(route_id) FROM routes",
                                    -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                sqlite3_close (db);
                return;
            }
            _routes.reserve (sqlite3_column_int (stmt, 0));
            sqlite3_finalize (stmt);

            if (sqlite3_prepare_v2 (db, "SELECT route_id FROM routes",
                                    -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                sqlite3_close (db);
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
            if (sqlite3_prepare_v2 (db, "SELECT count(trip_id) FROM trips",
                                    -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                sqlite3_close (db);
                return;
            }
            _trips.reserve (sqlite3_column_int (stmt, 0));
            sqlite3_finalize (stmt);

            if (sqlite3_prepare_v2 (db, "SELECT trip_id FROM trips",
                                    -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                sqlite3_close (db);
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
            if (sqlite3_prepare_v2 (db, "SELECT count(distinct shape_id) FROM shapes",
                                    -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                sqlite3_close (db);
                return;
            }
            _shapes.reserve (sqlite3_column_int (stmt, 0));
            sqlite3_finalize (stmt);

            if (sqlite3_prepare_v2 (db, "SELECT distinct shape_id FROM shapes",
                                    -1, &stmt, 0) != SQLITE_OK)
            {
                Rcpp::Rcerr << " x Can't prepare query\n  "
                    << sqlite3_errmsg (db) << "\n";
                sqlite3_finalize (stmt);
                sqlite3_close (db);
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

        sqlite3_close (db);


        
        //  --- testing
        Rcpp::Rcout << "\n *** test creating some things ...\n";

        // std::string tid ("1141101952-20180702170310_v67.28");
        
        // _agencies.emplace (std::piecewise_construct,
        //                    std::forward_as_tuple (aid),
        //                    std::forward_as_tuple (aid, this));
        // _routes.emplace (std::piecewise_construct,
        //                  std::forward_as_tuple (rid),
        //                  std::forward_as_tuple (rid, this));
        // _trips.emplace (std::piecewise_construct,
        //                 std::forward_as_tuple (tid),
        //                 std::forward_as_tuple (tid, this));


        Rcpp::Rcout << "\n\n *** Request information -> load ...\n";
        Rcpp::Rcout << "\n\n > Trip is run by "
            << _trips.begin ()->second.route ()
                    ->agency ()->agency_name ()
            << " and has "
            << _trips.begin ()->second.shape ()->path ().size ()
            << " points in its path\n";
    }

    std::string& Gtfs::dbname () 
    {
        return _dbname;
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


    /***************************************************** Agency */
    Agency::Agency (std::string& id, Gtfs* gtfs) : 
        gtfs (gtfs), _agency_id (id)
    {
        // Rcpp::Rcout << " + Create Agency " << _agency_id << "\n";
    }

    void Agency::load ()
    {
        sqlite3* db;
        if (sqlite3_open (gtfs->dbname ().c_str (), &db))
        {
            Rcpp::Rcerr << " x Unable to connect to database\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_close (db);
            return;
        }

        sqlite3_stmt* stmt;
        if (sqlite3_prepare_v2(db, "SELECT agency_name, agency_url, agency_phone, agency_timezone, agency_lang FROM agency WHERE agency_id=?",
                               -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            sqlite3_close (db);
            return;
        }
        if (sqlite3_bind_text (stmt, 1, _agency_id.c_str (),
                               -1, SQLITE_STATIC) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind agency id to query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            sqlite3_close (db);
            return; 
        }
        if (sqlite3_step (stmt) != SQLITE_ROW)
        {
            Rcpp::Rcerr << " x Couldn't get agency from db\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            sqlite3_close (db);
            return;
        }

        _agency_name = (char*)sqlite3_column_text (stmt, 0);
        _agency_url = (char*)sqlite3_column_text (stmt, 1);
        _agency_phone = (char*)sqlite3_column_text (stmt, 2);
        _agency_timezone = (char*)sqlite3_column_text (stmt, 3);
        _agency_lang = (char*)sqlite3_column_text (stmt, 4);

        sqlite3_finalize (stmt);
        sqlite3_close (db);

        loaded = true;
        Rcpp::Rcout << " + Agency " << _agency_id << " is loaded"
            << "\n   - Name: " << _agency_name 
            << "\n   - URL: " << _agency_url
            << "\n   - Phone: " << _agency_phone
            << "\n   - TZ: " << _agency_timezone
            << "\n   - Language: " << _agency_lang << "\n";

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
        if (!loaded) load();
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
        sqlite3* db;
        if (sqlite3_open (gtfs->dbname ().c_str (), &db))
        {
            Rcpp::Rcerr << " x Unable to connect to database\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_close (db);
            return;
        }

        sqlite3_stmt* stmt;
        if (sqlite3_prepare_v2(db, "SELECT route_short_name, route_long_name, route_type, agency_id, version FROM routes WHERE route_id=?",
                               -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            sqlite3_close (db);
            return;
        }
        if (sqlite3_bind_text (stmt, 1, _route_id.c_str (),
                               -1, SQLITE_STATIC) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind route id to query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            sqlite3_close (db);
            return; 
        }
        if (sqlite3_step (stmt) != SQLITE_ROW)
        {
            Rcpp::Rcerr << " x Couldn't get route from db\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            sqlite3_close (db);
            return;
        }

        _route_short_name = (char*)sqlite3_column_text (stmt, 0);
        _route_long_name = (char*)sqlite3_column_text (stmt, 1);
        _route_type = sqlite3_column_int (stmt, 2);
        std::string agency_id = (char*)sqlite3_column_text (stmt, 3);
        _version = (float)sqlite3_column_double (stmt, 4);

        _agency = gtfs->find_agency (agency_id);

        sqlite3_finalize (stmt);
        sqlite3_close (db);

        loaded = true;
        Rcpp::Rcout << " + Route " << _route_id << " is loaded"
            << "\n   - Number: " << _route_short_name
            << " " << _route_long_name
            << "\n   - Agency: " 
            << (_agency != nullptr ? _agency->agency_name () : "null")
            << "\n   - Version: " << _version << "\n";
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
        if (!loaded) load();
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
        // Rcpp::Rcout << " + Create Trip " << _trip_id << "\n";
    }

    void Trip::load ()
    {
        sqlite3* db;
        if (sqlite3_open (gtfs->dbname ().c_str (), &db))
        {
            Rcpp::Rcerr << " x Unable to connect to database\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_close (db);
            return;
        }

        sqlite3_stmt* stmt;
        if (sqlite3_prepare_v2(db, "SELECT route_id, shape_id, service_id, block_id, direction_id, trip_headsign, version FROM trips WHERE trip_id=?",
                               -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            sqlite3_close (db);
            return;
        }
        if (sqlite3_bind_text (stmt, 1, _trip_id.c_str (),
                               -1, SQLITE_STATIC) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind trip id to query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            sqlite3_close (db);
            return; 
        }
        if (sqlite3_step (stmt) != SQLITE_ROW)
        {
            Rcpp::Rcerr << " x Couldn't get trip from db\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            sqlite3_close (db);
            return;
        }

        std::string routeid = (char*)sqlite3_column_text (stmt, 0);
        _route = gtfs->find_route (routeid);
        std::string shapeid = (char*)sqlite3_column_text (stmt, 1);
        _shape = gtfs->find_shape (shapeid);
        std::string serviceid = (char*)sqlite3_column_text (stmt, 2);
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
        sqlite3_close (db);

        loaded = true;
        Rcpp::Rcout << " + Trip " << _trip_id << " is loaded"
            << "\n   - Route: "
            << (_route != nullptr ? _route->route_short_name () : "missing")
            << "\n   - Headsign: " << _trip_headsign
            << "\n   - Version: " << _version << "\n";
    }

    void Trip::unload () { unload (false); }

    void Trip::unload (bool complete)
    {
        completed = complete;
        loaded = false;
        _route = nullptr;
        _shape = nullptr;
        // _calendar = nullptr;
        _block_id = "";
        _trip_headsign = "";
        Rcpp::Rcout << " + Trip " << _trip_id << " is unloaded\n";
    }

    std::string& Trip::trip_id () { 
        if (!loaded) load ();
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
    // Calendar* Trip::calendar () { 
    //     if (!loaded) load ();
    //     return &_calendar; 
    // }
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
        sqlite3* db;
        if (sqlite3_open (gtfs->dbname ().c_str (), &db))
        {
            Rcpp::Rcerr << " x Unable to connect to database\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_close (db);
            return;
        }

        sqlite3_stmt* stmt;
        if (sqlite3_prepare_v2 (db, "SELECT count(shape_id) FROM shapes WHERE shape_id=?",
                                -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            sqlite3_close (db);
            return;
        }
        if (sqlite3_bind_text (stmt, 1, _shape_id.c_str (),
                               -1, SQLITE_STATIC) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind trip id to query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            sqlite3_close (db);
            return; 
        }
        if (sqlite3_step (stmt) != SQLITE_ROW)
        {
            Rcpp::Rcerr << " x Couldn't get row count from db\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            sqlite3_close (db);
            return;
        }
        _path.reserve (sqlite3_column_int (stmt, 0));
        sqlite3_finalize (stmt);

        if (sqlite3_prepare_v2(db, "SELECT shape_pt_lat, shape_pt_lon, shape_dist_traveled, version FROM shapes WHERE shape_id=? ORDER BY shape_pt_sequence",
                               -1, &stmt, 0) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't prepare query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            sqlite3_close (db);
            return;
        }
        if (sqlite3_bind_text (stmt, 1, _shape_id.c_str (),
                               -1, SQLITE_STATIC) != SQLITE_OK)
        {
            Rcpp::Rcerr << " x Can't bind shape id to query\n  "
                << sqlite3_errmsg (db) << "\n";
            sqlite3_finalize (stmt);
            sqlite3_close (db);
            return; 
        }
        
        while (sqlite3_step (stmt) == SQLITE_ROW)
        {
            _path.emplace_back (sqlite3_column_double (stmt, 0),
                                sqlite3_column_double (stmt, 1),
                                sqlite3_column_double (stmt, 2));
        }
        _version = (float)sqlite3_column_double (stmt, 3);
        
        sqlite3_finalize (stmt);
        sqlite3_close (db);

        loaded = true;
        Rcpp::Rcout << " + Shape " << _shape_id << " is loaded"
            << "\n   - Path: " << _path.size () << " coordinates"
            << "\n   - Version: " << _version << "\n";
    }

    void Shape::unload () { unload (false); }

    void Shape::unload (bool complete)
    {
        completed = complete;
        loaded = false;
        Rcpp::Rcout << " + Shape " << _shape_id << " is unloaded\n";
    }

    std::string& Shape::shape_id () { 
        if (!loaded) load ();
        return _shape_id; 
    }

    std::vector<ShapePt>& Shape::path ()
    {
        if (!loaded) load ();
        return _path;
    }

    float Shape::version () { 
        if (!loaded) load ();
        return _version; 
    }


}; // namespace Gtfs