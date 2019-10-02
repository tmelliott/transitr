/**
 * Define the ETA methods.
 */

#include "etas.h"
#if VERBOSE == 2
#include "timing.h"
#endif


namespace Gtfs {

    void Trip::update (uint64_t& t, RNG& rng)
    {
        // Update trip state
        if (!loaded) load ();
    }

    etavector Trip::get_etas ()
    {
        etavector etas;
        return etas;
    }

} // namespace Gtfs
