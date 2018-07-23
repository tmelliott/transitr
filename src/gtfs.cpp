#include "gtfs.h"
#include "vendor/sqlite3/sqlite3.h"

namespace Gtfs 
{

    /***************************************************** GTFS */
    Gtfs::Gtfs (std::string& name) : _dbname (name) 
    {
        Rcpp::Rcout << "Connected to GTFS database `"
            << _dbname << "`\n";

        // Now load all of the things ...
        Rcpp::Rcout << "\n *** test creating some things ...\n";

        std::string aid ("HE");
        std::string rid ("07005-20180702170310_v67.28");
        std::string tid ("1141101952-20180702170310_v67.28");
        
        _agencies.emplace (std::piecewise_construct,
                           std::forward_as_tuple (aid),
                           std::forward_as_tuple (aid, this));
        _routes.emplace (std::piecewise_construct,
                         std::forward_as_tuple (rid),
                         std::forward_as_tuple (rid, this));
        _trips.emplace (std::piecewise_construct,
                        std::forward_as_tuple (tid),
                        std::forward_as_tuple (tid, this));


        Rcpp::Rcout << "\n\n *** test loading some things ...\n";
        // _agencies.begin ()->second.load ();
        // _routes.begin ()->second.load ();
        _trips.begin ()->second.load ();

        Rcpp::Rcout << "\n\n *** test un-loading some things ...\n";
        _agencies.begin ()->second.unload ();
        _routes.begin ()->second.unload ();
        _trips.begin ()->second.unload ();
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


    /***************************************************** Agency */
    Agency::Agency (std::string& id, Gtfs* gtfs) : 
        _agency_id (id), gtfs (gtfs)
    {
        Rcpp::Rcout << " + Create Agency " << _agency_id << "\n";
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
        _route_id (id), gtfs (gtfs)
    {
        Rcpp::Rcout << " + Create Route " << _route_id << "\n";
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
        completed = complete;
        loaded = false;
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
        _trip_id (id), gtfs (gtfs)
    {
        Rcpp::Rcout << " + Create Trip " << _trip_id << "\n";
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
        Rcpp::Rcout << " + Trip " << _trip_id << " is unloaded\n";
    }

    std::string& Trip::trip_id () { return _trip_id; }
    Route* Trip::route () { return _route; }
    // Shape* Trip::shape () { return &_shape; }
    // Calendar* Trip::calendar () { return &_calendar; }
    std::string& Trip::block_id () { return _block_id; }
    bool Trip::direction_id () { return _direction_id; }
    std::string& Trip::trip_headsign () { return _trip_headsign; }
    float Trip::version () { return _version; }


}; // namespace Gtfs