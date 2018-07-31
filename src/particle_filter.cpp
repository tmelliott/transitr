/**
 * Define the particle filter functions.
 */
#include "particle_filter.h"

namespace Gtfs {

    void Vehicle::initialize ()
    {
        _state.reserve (_N);
        double dmax = _trip->shape ()->path ().back ().distance;
        for (int i=0; i<_N; ++i)
        {
            // initialize each particle with a state X(i)
            
        }
    }

    void Vehicle::mutate ()
    {

    }

}; // namespace Gtfs
