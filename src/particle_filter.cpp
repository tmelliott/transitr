/**
 * Define the particle filter functions.
 */
#include "particle_filter.h"

namespace Gtfs {

    void Vehicle::initialize ()
    {
        _state.reserve (_N);
        double dmax = _trip->shape ()->path ().back ().distance;
        // std::cout << "\n + Initialize " << _vehicle_id << ": ";
        for (int i=0; i<_N; ++i)
        {
            // initialize each particle with a state X(i)
            _state.emplace_back (10.0, 1.0, this);
            // std::cout << " [" << _state.back ().get_distance () << "]";
        }
    }

    void Vehicle::mutate ()
    {

    }



    Particle::Particle (double d, double s, Vehicle* v)
    {
        vehicle = v;
        distance = d;
        speed = s;
    }

    double Particle::get_distance ()
    {
        return distance;
    }

}; // namespace Gtfs
