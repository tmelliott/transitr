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

            std::cout << "\n\n * Summarising ..." << std::endl;
            arrival_times = get_etas ();

            print_etas ();
        }

        std::cout << "\n\n------------------------------------------\n";

    }

    /**
     * Forecast arrival times using normal theory
     * @param rng an RNG
     */
    void Trip::forecast (RNG& rng)
    {
        Eigen::IOFormat intMat (
            Eigen::StreamPrecision, 0, ", ", "\n", "  [", "]"
        );

        int M = _stops.size ();
        auto segs = _shape->segments ();
        int L = segs.size ();
        Eigen::MatrixXi Hseg (Eigen::MatrixXi::Zero (L, M-1));
        int m=0;
        for (int l=0; l<L; l++)
        {
            if (_stops.at (m+1).distance <= segs.at (l).distance) m++;
            Hseg (l, m) = 1;
        }
        // std::cout << "\n a big matrix ...\n" << Hseg.format (intMat) << "\n";


        if (gtfs->parameters ()->eta_model == 0 && _vehicle != nullptr)
        {
            // iterate over vehicle state
            int N (_vehicle->state ()->size ());
            N = std::min (N, 10);
            Eigen::MatrixXi seg_tt (Eigen::MatrixXi::Zero (N, L));
            Eigen::MatrixXi dwell_t (Eigen::MatrixXi::Zero (N, M-1));

            Particle* p;
            int l;
            int tt;
            double dt;
            double v;
            //  vehicle speed too slow?
            bool use_particle_speed (_vehicle->speed () < 2);
            for (int i=0; i<N; i++)
            {
                tt = 0;
                p = &(_vehicle->state ()->at (i));
                // time to end of current segment
                l = find_segment_index (p->get_distance (), &segs);
                v = use_particle_speed ? p->get_speed () : 
                    segs.at (l).segment->sample_speed (rng);
                tt += 
                    round ((
                        segs.at (l).distance + 
                        segs.at (l).segment->length () - 
                        p->get_distance ()
                    ) / p->get_speed ());

                seg_tt (i, l) = fmin (tt, 60*30); // no more than 30 mins!

                l++;
                while (l < L)
                {
                    // fetch tt from scheduled travel time ...
                    // tt += _stops.at (l+1).arrival_time - _stops.at (l).departure_time;
                    tt += segs.at (l).segment->sample_travel_time (rng, tt);
                    seg_tt (i, l) = fmin (tt, 60*60*2); // no more than 2 hours!!
                    l++;
                }

                int m = find_stop_index (p->get_distance (), &(_stops));
                for (m; m<M-1; m++)
                {
                    // first and last column are zeros
                    if (m == 0) continue;
                    // if bus has arrived but not departed, (and GPS is at stop),
                    // then include dwell time, otherwise continue ()
                    if (rng.runif () < gtfs->parameters ()->pr_stop)
                    {
                        dt = -1.0;
                        while (dt <= 0 || dt > 5*60)
                        {
                            dt = rng.rnorm () * gtfs->parameters ()->dwell_time_var +
                                gtfs->parameters ()->dwell_time;
                        }
                        dwell_t (i, m) = round (gtfs->parameters ()->gamma + dt);
                        if (i > 0) dwell_t (i, m) += dwell_t (i-1, m);
                    }
                }
            }
            std::cout << "\n" << seg_tt.format (intMat) << "\n";

            // Convert segment mat to link mat
            Eigen::MatrixXi link_tt = seg_tt * Hseg;

            std::cout << "\n" << link_tt.format (intMat) << "\n";
            std::cout << "\n" << dwell_t.format (intMat) << "\n";

            _eta_matrix = link_tt + dwell_t;
            std::cout << "\n" << _eta_matrix.format (intMat) << "\n";
        }
        else
        {

        }
    }

    etavector Trip::get_etas ()
    {
        etavector etas;
        int M (_stops.size ());
        etas.resize (M);
        int N (_eta_matrix.rows ());

        Eigen::VectorXi col_m;
        for (int m=_stop_index;m<M-1;m++)
        {
            col_m = _eta_matrix.col (m);
            if (col_m.isZero ()) continue;

            std::sort (col_m.data (), col_m.data () + col_m.size ());
            // median
            etas.at (m+1).stop_id = _stops.at (m).stop->stop_id ();
            etas.at (m+1).estimate = _timestamp + col_m (floor (N/2));
        }

        return etas;
    }

    void Trip::print_etas ()
    {
        // if (!state_initialised) return;

        std::cout << "\n\n - ETAs for route "
            << route ()->route_short_name ()
            << " ("
            << stops ().at (0).arrival_time
            << ") with " << stops ().size () << " stops\n"
            << "              schedule       eta  delay  (error)";

        // std::string
        for (int i=0; i<stops ().size (); i++)
        {
            if (arrival_times.at (i).estimate == 0) continue;

            std::cout << "\n   + stop "
                << std::setw (2) << (i+1) << ": "
                << stops ().at (i).arrival_time << "  "
                << Time (arrival_times.at (i).estimate) << "  "
                << std::setw (5)
                << (arrival_times.at (i).estimate - stops ().at (i).arrival_time)
                << "  ("
                // << uncertainty.at (i)
                << ")";
        }
    }

} // namespace Gtfs
