#include "realtime_feed.h"
#include <regex>

#include "timing.h"

static size_t WriteCallback(void *contents, size_t size, size_t nmemb, void *userp)
{
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}

RealtimeFeed::RealtimeFeed (std::vector<std::string>& urls, List& hdrs) : _urls (urls)
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

/**
 * Update the protobuf realtime feed with the provided URL(s)
 *
 * Each URL will be iteratively added to the feed - trip updates
 * and vehicle positions will be combined into a single feed.
 * 
 * @return integer return code: 0 - ok, 1 - curl failed, 2 - parse feed failed, 
 *                 5 - simuation has completed
 */
int RealtimeFeed::update ()
{
    _feed.Clear ();
    _n_vehicles = 0;
    _n_trip_updates = 0;

    for (auto url = _urls.begin (); url != _urls.end (); ++url)
    {
        CURL *curl;
        CURLcode res;
        curl_global_init (CURL_GLOBAL_DEFAULT);
        curl = curl_easy_init ();

        std::string readBuffer;
        if (curl)
        {
            // add the headers
            struct curl_slist *chunk = NULL;

            for (auto header: _headers)
            {
                const char* hstr = header.c_str ();
                chunk = curl_slist_append (chunk, hstr);
            }

            res = curl_easy_setopt(curl, CURLOPT_HTTPHEADER, chunk);
            // curl_easy_setopt(curl, CURLOPT_VERBOSE, 1L);

            const char* urlstr = url->c_str ();
            curl_easy_setopt (curl, CURLOPT_URL, urlstr);
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

        if (std::regex_match (readBuffer, std::regex ("(simulation complete)"))) {
            return 5;
        }

        std::istringstream buf (readBuffer);
        transit_realtime::FeedMessage feed;
        if (!feed.ParseFromIstream (&buf)) {
            Rcerr << "\n x Failed to parse GTFS realtime feed!\n";
            return 2;
        }

        if (!_feed.has_header ())
        {
            // need to add the header
            transit_realtime::FeedHeader* header;
            header = _feed.mutable_header ();
            header->set_gtfs_realtime_version ("2.0");
            header->set_timestamp (feed.header ().timestamp ());
        }

        // copy entities into _feed
        for (int i=0; i<feed.entity_size (); ++i)
        {
            // auto ent = feed->entity (i);
            // auto ent = _feed.add_entity ();
            // ent = feed.entity (i);
            _feed.add_entity ()->CopyFrom (feed.entity (i));
            if (feed.entity (i).has_vehicle ()) _n_vehicles++;
            if (feed.entity (i).has_trip_update ()) _n_trip_updates++;
        }
    }

    return 0;
}

transit_realtime::FeedMessage* RealtimeFeed::feed ()
{
    return &_feed;
}

void load_vehicles (Gtfs::vehicle_map* vehicles,
                    transit_realtime::FeedMessage* feed,
                    Gtfs::Gtfs* gtfs, Gtfs::par* params)
{

    for (int i=0; i<feed->entity_size (); ++i)
    {
        auto ent = feed->entity (i);
        if (ent.has_vehicle () && ent.vehicle ().has_vehicle ())
        {
            // std::cout << " - vehicle " << ent.vehicle ().vehicle ().id ()
                // << " - " << ent.vehicle ().trip ().id () << "\n";
            std::string id (ent.vehicle ().vehicle ().id ());
            // if (id != "2C9F") continue;
            auto vs = vehicles->find (id);
            if (vs == vehicles->end ())
            {
                auto r = vehicles->emplace (std::piecewise_construct,
                                            std::forward_as_tuple (id), 
                                            std::forward_as_tuple (id, params));
                if (r.second)
                {
                    r.first->second.update (ent.vehicle (), gtfs);
                }
            }
            else
            {
                vs->second.update (ent.vehicle (), &(*gtfs));
            }
        } // end if vehicle position
        
        if (ent.has_trip_update () && ent.trip_update ().stop_time_update ().size () > 0)
        {
            std::string id (ent.trip_update ().vehicle ().id ());
            // if (id != "2C9F") continue;
            auto vs = vehicles->find (id);
            if (vs == vehicles->end ())
            {
                auto r = vehicles->emplace (std::piecewise_construct,
                                            std::forward_as_tuple (id), 
                                            std::forward_as_tuple (id, params));
                if (r.second)
                {
                    r.first->second.update (ent.trip_update (), gtfs);
                }
            }
            else
            {
                vs->second.update (ent.trip_update (), &(*gtfs));
            }
        } // end if trip update
    } // end feed->entity

    // update all the vehicles
    for (auto v = vehicles->begin (); v != vehicles->end (); ++v)
    {
        v->second.update (gtfs);
    }

}

void write_vehicles (Gtfs::vehicle_map* vehicles, std::string& file)
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
        if (v->second.complete () || !v->second.valid ()) continue;
        transit_realtime::FeedEntity* entity = feed.add_entity ();
        entity->set_id (v->second.vehicle_id ());

        // --- Vehicle Position
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

        // --- Trip Update (ETAs)
        transit_realtime::TripUpdate* tu = entity->mutable_trip_update ();
        if (v->second.trip () != nullptr)
        {
            transit_realtime::TripDescriptor* trip = tu->mutable_trip ();
            trip->set_trip_id (v->second.trip ()->trip_id ());
            if (v->second.trip ()->route () != nullptr) 
            {
                trip->set_route_id (v->second.trip ()->route ()->route_id ());
            }
        }
        vehicle = tu->mutable_vehicle ();
        vehicle->set_id (v->second.vehicle_id ());

        // Stop Time Events
        Gtfs::etavector etas (v->second.get_etas ());
        for (int si=0; si<etas.size (); ++si)
        {
            transit_realtime::TripUpdate::StopTimeUpdate* stu = tu->add_stop_time_update ();
            stu->set_stop_sequence (si+1);
            transit_network::TimePrediction* tpi = stu->MutableExtension(transit_network::eta);
            tpi->set_estimate (etas.at (si).estimate);
            for (auto q : etas.at (si).quantiles)
            {
                transit_network::Quantile* qi = tpi->add_quantiles ();
                qi->set_quantile (q.quantile);
                qi->set_value (q.time);
            }
        }

    }

    // write the feed to a file
    std::fstream output (file.c_str (), std::ios::out | std::ios::trunc | std::ios::binary);
    if (!feed.SerializeToOstream (&output)) 
    {
        std::cerr << "\n x Failed to write feed to `" << file << "`\n";
    }

}
