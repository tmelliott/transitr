#include "realtime_feed.h"

#include "timing.h"

static size_t WriteCallback(void *contents, size_t size, size_t nmemb, void *userp)
{
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}

RealtimeFeed::RealtimeFeed (std::string& url, List& hdrs) : _url (url)
{
    _headers.reserve (hdrs.size ());
    for (int i=0; i<hdrs.size (); i++)
    {
        List hdr = hdrs[i];
        String hdrn = hdr["name"];
        String hdrv = hdr["value"];
        std::ostringstream hss;
        hss << (std::string)hdrn << ": " << (std::string)hdrv;
        _headers.push_back (hss.str ());
    }
}

int RealtimeFeed::update ()
{
    CURL *curl;
    CURLcode res;
    curl_global_init (CURL_GLOBAL_DEFAULT);
    curl = curl_easy_init ();

    _feed.Clear ();

    std::string readBuffer;
    if (curl)
    {
        // add the headers
        struct curl_slist *chunk = NULL;

        for (auto header: _headers)
        {
            chunk = curl_slist_append (chunk, header.c_str ());
        }

        res = curl_easy_setopt(curl, CURLOPT_HTTPHEADER, chunk);
        // curl_easy_setopt(curl, CURLOPT_VERBOSE, 1L);

        curl_easy_setopt (curl, CURLOPT_URL, _url.c_str ());
        curl_easy_setopt (curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt (curl, CURLOPT_WRITEDATA, &readBuffer);
        res = curl_easy_perform (curl);
        if (res != CURLE_OK)
        {
            Rcerr << "curl_easy_perform() failed\n";
            return 1;
            // curl_easy_strerror(res));
        }
        curl_easy_cleanup (curl);
        curl_slist_free_all (chunk);
    }
    curl_global_cleanup ();

    std::istringstream buf (readBuffer);
    if (!_feed.ParseFromIstream (&buf)) {
        Rcerr << "\n x Failed to parse GTFS realtime feed!\n";
        return 2;
    }

    return 0;
}

transit_realtime::FeedMessage* RealtimeFeed::feed ()
{
    return &_feed;
}

void load_vehicles (Gtfs::vehicle_map* vehicles,
                    transit_realtime::FeedMessage* feed,
                    Gtfs::Gtfs* gtfs, int n, double err)
{
#if VERBOSE == 2
    Timer timer;
#endif
    for (int i=0; i<feed->entity_size (); ++i)
    {
        auto ent = feed->entity (i);
        if (!ent.has_vehicle ()) continue;
        if (!ent.vehicle ().has_vehicle ()) continue;
        
        std::string id (ent.vehicle ().vehicle ().id ());
#if VERBOSE == 2
        std::cout << " + loading vehicle " << id;
#endif
        auto vs = vehicles->find (id);
#if VERBOSE == 2
        std::cout << " (" << timer.cpu_seconds () << "ms)";
        timer.reset ();
#endif
        if (vs == vehicles->end ())
        {
#if VERBOSE == 2
            std::cout << " - insert";
#endif
            auto r = vehicles->emplace (std::piecewise_construct,
                                        std::forward_as_tuple (id), 
                                        std::forward_as_tuple (id, n, err));
#if VERBOSE == 2
            std::cout << " (" << timer.cpu_seconds () << "ms)";
            timer.reset ();
#endif
            if (r.second)
            {
                r.first->second.update (ent.vehicle (), gtfs);
            }
        }
        else
        {
            vs->second.update (ent.vehicle (), &(*gtfs));
        }
#if VERBOSE == 2
        std::cout << " => TOTAL UPDATE (" << timer.cpu_seconds () << "ms)\n";
        timer.reset ();
#endif

        // if (i >= 10) break;
    }
}

void write_vehicles (Gtfs::vehicle_map* vehicles)
{
    // create a new feed
    transit_realtime::FeedMessage feed;

    // set the header information
    transit_realtime::FeedHeader* header;
    header = feed.mutable_header ();
    header->set_gtfs_realtime_version ("2.0");
    std::time_t curtime = std::time (nullptr);
    header->set_timestamp (curtime);

    // write vehicles
    for (auto v = vehicles->begin (); v != vehicles->end (); ++v)
    {
        if (!v->second.valid ()) continue;
        transit_realtime::FeedEntity* entity = feed.add_entity ();
        entity->set_id (v->second.vehicle_id ());

        transit_realtime::VehiclePosition* vp = entity->mutable_vehicle ();
        if (v->second.trip () != nullptr) 
        {
            transit_realtime::TripDescriptor* trip = vp->mutable_trip ();
            trip->set_trip_id (v->second.trip ()->trip_id ());
            if (v->second.trip ()->route () != nullptr) 
            {
                trip->set_route_id (v->second.trip ()->route ()->route_id ());
            }
        }
        transit_realtime::VehicleDescriptor* vehicle = vp->mutable_vehicle ();
        vehicle->set_id (v->second.vehicle_id ());

        // the observed position
        transit_realtime::Position* position = vp->mutable_position ();
        position->set_latitude (v->second.position ().latitude);
        position->set_longitude (v->second.position ().longitude);

        // the modeled position
        transit_realtime::Position* 
        position_estimate = vp->MutableExtension (transit_network::position_estimate);
        double dbar = v->second.distance ();
        if (v->second.trip ()->shape () != nullptr)
        {
            latlng vpos = v->second.trip ()->shape ()->coordinates_of (dbar);
            position_estimate->set_latitude (vpos.latitude);
            position_estimate->set_longitude (vpos.longitude);
        }
        position_estimate->set_odometer (dbar);
        position_estimate->set_speed (v->second.speed ());

    }

    // write the feed to a file
    std::fstream output ("at_predictions.pb", 
                         std::ios::out | std::ios::trunc | std::ios::binary);
    if (!feed.SerializeToOstream (&output)) {
        std::cerr << "\n x Failed to write feed.\n";
    }

}
