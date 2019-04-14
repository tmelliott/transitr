#include "gtfs.h"

#if SIMULATION
namespace Gtfs
{
    void Vehicle::store_state (std::string type)
    {
        if (!store_created)
        {
            f.open (store_name.c_str ());
            f
                << "vehicle_id"
                << ",trip_id"
                << ",state_type"
                << ",event_timestamp"
                << ",event_type"
                << ",event_latitude"
                << ",event_longitude"
                << ",event_dist_to_route"
                << ",event_dist_along_route"
                << ",event_stop_index"
                << ",vehicle_timestamp"
                << ",delta"
                << ",vehicle_distance"
                << ",vehicle_speed"
                << ",vehicle_latitude"
                << ",vehicle_longitude"
                << ",vehicle_position_error"
                << ",sum_llh"
                << ",Neff"
                << ",action"
                << "\n";
            f.close ();
            store_created = true;
        }

        f.open (store_name.c_str (), std::ofstream::app);

        // current "event"
        auto& e = time_events.at (current_event_index);
        double d = _trip->shape ()->distance_of (e.position);
        latlng px = _trip->shape ()->coordinates_of (d);

        double llh = std::accumulate (
            _state.begin (), 
            _state.end (), 
            0.0,
            [](double a, Particle& p) 
            {
                return a + p.get_weight () * exp (p.get_ll ());
            }
        );

        f
            << _vehicle_id
            << "," << _trip->trip_id ()
            << "," << type
            << "," << e.timestamp
            << "," << e.type_name ()
            << "," << (e.type == EventType::gps ? std::to_string (e.position.latitude) : "")
            << "," << (e.type == EventType::gps ? std::to_string (e.position.longitude) : "")
            << "," << (e.type == EventType::gps ? std::to_string (distanceEarth (e.position, px)) : "")
            << "," << (e.type == EventType::gps ? std::to_string (_trip->shape ()->distance_of (e.position)) : "")
            << "," << (e.type == EventType::gps ? "" : std::to_string (e.stop_index))
            << "," << _timestamp
            << "," << _delta
            << "," << this->distance ()
            << "," << this->speed ()
            << "," << _position.latitude
            << "," << _position.longitude
            << "," << (e.type == EventType::gps ? std::to_string (distanceEarth (e.position, _position)) : "")
            << "," << llh
            << "," << _Neff
            << "," << action
            << "\n";

        f.close ();

        action = "";
    }
}
#endif
