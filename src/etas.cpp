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

        std::cout << "\n * Updating trip " << _trip_id;

        int eta_model (gtfs->parameters ()->eta_model);
        std::cout
            << "\n   - Estimating ETAs using Model " << eta_model;

        std::cout << "\n   - Vehicle: "
            << (_vehicle == nullptr ? "none" : _vehicle->vehicle_id ());


        std::cout << "\n\n------------------------------------------\n";

    }

    etavector Trip::get_etas ()
    {
        etavector etas;
        return etas;
    }

} // namespace Gtfs
