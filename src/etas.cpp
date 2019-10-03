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

        std::cout << "\n * Updating trip " << _trip_id
            << " (" << _route->route_short_name () << ")";

        int eta_model (gtfs->parameters ()->eta_model);
        std::cout
            << "\n   - Estimating ETAs using Model " << eta_model;

        std::cout << "\n   - Vehicle: "
            << (_vehicle == nullptr ? "none" : _vehicle->vehicle_id ());

        if (_vehicle != nullptr)
        {
            // has vehicle been updated?
            if (_vehicle->timestamp () > _timestamp)
            {
                std::cout << " (updated - ";

                Event* event = _vehicle->latest_event ();

                _segment_progress = 0;
                if (event->type == EventType::gps)
                {
                    std::cout << "gps";
                    _event_type = 0;
                    double d = _vehicle->distance ();
                    _segment_index = find_segment_index (d, &(_shape->segments ()));
                    _stop_index = find_stop_index (d, &(_stops));
                    double Dl = _shape->segments ().at (_segment_index).distance;
                    double Ll = _shape->segments ().at (_segment_index).segment->length ();
                    _segment_progress = (d - Dl) / Ll;
                } 
                else
                {
                    if (event->type == EventType::arrival)
                    {
                        _event_type = 1;
                        std::cout << "arrival";
                    }
                    else
                    {
                        _event_type = 2;
                        std::cout << "departure";
                    }

                    _stop_index = event->stop_index;
                    double d = _stops.at (_stop_index).distance;
                    _segment_index = find_segment_index (d, &(_shape->segments ()));
                }
                std::cout << ")";
                
            }
        }
        _timestamp = t;

        std::cout
            << "\n   - Event Type: " << _event_type
            << "\n   - Stop index: " << (_stop_index+1)
            << " of " << _stops.size ()
            << "\n   - Segment index: " << (_segment_index+1)
            << " of " << _shape->segments ().size ()
            << " (" << (_segment_progress*100) << "%)";

        if (_stop_index < _shape->segments ().size ())
        {
            std::cout << "\n\n * Forecasting ..." << std::endl;
            forecast (rng);
        }

        std::cout << "\n\n------------------------------------------\n";

    }

    /**
     * Forecast arrival times using normal theory
     * @param rng an RNG
     */
    void Trip::forecast (RNG& rng)
    {
        Eigen::IOFormat intMat (Eigen::StreamPrecision, 0, ", ", "\n", "  [", "]");
        if (gtfs->parameters ()->eta_model == 0 && _vehicle != nullptr)
        {
            // iterate over vehicle state
            int N (_vehicle->state ()->size ());
            N = std::min (N, 10);
            int L (_shape->segments ().size ());
            Eigen::MatrixXd seg_tt (Eigen::MatrixXd::Zero (N, L));

            auto segs = _shape->segments ();

            Particle* p;
            int l;
            int tt;
            for (int i=0; i<N; i++)
            {
                tt = 0;
                p = &(_vehicle->state ()->at (i));
                // time to end of current segment
                l = find_segment_index (p->get_distance (), &segs);
                tt += 
                    round ((
                        segs.at (l).distance + 
                        segs.at (l).segment->length () - 
                        p->get_distance ()
                    ) / p->get_speed ());
                seg_tt (i, l) = tt;

                l++;
                while (l < L)
                {
                    tt += segs.at (l).segment->sample_travel_time (rng, tt);
                    seg_tt (i, l) = tt;
                    l++;
                }
            }
            std::cout << "\n" << seg_tt.format (intMat) << "\n";
        }
        else
        {

        }
    }

    etavector Trip::get_etas ()
    {
        etavector etas;
        return etas;
    }

} // namespace Gtfs
