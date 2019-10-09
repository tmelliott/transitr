#include <string>
#include <iostream>
#include <fstream>

#include "gtfs.h"

using namespace Rcpp;

void processFile (std::string& f, std::fstream& o, Gtfs::Gtfs* gtfs)
{
    std::fstream input (f, std::ios::in | std::ios::binary);
    transit_realtime::FeedMessage feed;
    if (!feed.ParseFromIstream (&input)) {
        Rcerr << "\n x Failed to parse GTFS realtime file:\n" << f;
    }

    transit_realtime::FeedEntity* ent;
    transit_realtime::TripUpdate* tu;

    std::string tid;
    Gtfs::Trip* trip;
    uint64_t t;

    for (int i=0; i<feed.entity_size (); ++i)
    {
        ent = feed.mutable_entity (i);
        tu = ent->mutable_trip_update ();
        if (tu->timestamp() == 0) continue;

        tid = tu->trip ().trip_id ();
        trip = gtfs->find_trip (tid);
        t = tu->timestamp ();
        int delay;
        for (auto stu : tu->stop_time_update ())
        {
            if (stu.HasExtension (transit_network::current_delay))
            {
                delay = stu.GetExtension (transit_network::current_delay);
            }
            else if (stu.has_arrival ())
            {
                delay = stu.arrival ().delay ();
            } 
            else if (stu.has_departure ())
            {
                delay = stu.departure ().delay ();
            }
            else
            {
                delay = 0;
            }
            o 
                << tu->trip ().trip_id () 
                << "," << tu->trip ().route_id ()
                << "," << (tu->has_vehicle () ? tu->vehicle ().id () : "")
                << "," << t
                << "," << stu.stop_sequence ()
                << "," << delay
                << "," << (stu.HasExtension (transit_network::eta) ? stu.GetExtension (transit_network::eta).estimate () : 0)
                << "," << trip->stops ().at (stu.stop_sequence ()-1).arrival_time.asUNIX (t)
                << "\n";
        }
    }
}

// [[Rcpp::export]]
void processEtas (StringVector files, StringVector out, StringVector dbname)
{
    GOOGLE_PROTOBUF_VERIFY_VERSION;

    std::string db = as<std::string> (dbname[0]);
    Gtfs::Gtfs gtfs (db);

    std::string file;
    std::string outfile (out[0]);
    std::fstream o;
    o.open (outfile, std::ios::out);
    int i=0;
    int N=files.size ();
    int perc, newperc;
    std::cout << "\n";
    for (auto f=files.begin (); f != files.end (); ++f)
    {
        newperc = (++i * 20 / N) * 5;
        if (newperc > perc)
        {
            perc = newperc;
            std::cout << "\r * Processing " << N << " files ... " << perc << "%";
        }
        file = as<std::string> (*f);
        processFile (file, o, &gtfs);
    }
    o.close ();
    std::cout << " ... done!\n";

    // DataFrame result = 
    //     DataFrame::create(
    //         // _["trip_id"] = tripIds,
    //         // _["route_id"] = routeIds,
    //         // _["timestamp"] = dt,
    //         // _["stop_index"] = stopIndices,
    //         // _["arrival_time"] = arrivalTimes,
    //         // _["stringsAsFactors"] = false
    //     );

    // // return result as a tibble
    // result.attr ("class") = 
    //     StringVector::create ("tbl_df", "tbl", "data.frame");

}

