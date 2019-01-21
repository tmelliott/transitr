/**
 * Define the ETA methods.
 */

#include "etas.h"
#if VERBOSE == 2
#include "timing.h"
#endif


namespace Gtfs {


    void Vehicle::predict_etas (RNG& rng)
    {
        if (!valid ()) return;

#if VERBOSE == 2
        Timer timer;
        std::cout << "\n- vehicle " << _vehicle_id << " - predicting etas";
#endif
        for (auto p = _state.begin (); p != _state.end (); ++p) 
        {
            p->predict_etas (rng);
#if VERBOSE == 2
            std::cout << "\n       => ";
            for (auto eta : p->get_arrival_times ()) {
                if (eta <= _timestamp) std::cout << "*, ";
                else std::cout << ((eta - _timestamp) / 60) << ", ";
            }
#endif
            // std::cout << "\n VERSUS ... ";
            // for (int i=0; i<p->get_arrival_times ().size (); ++i) 
            //     std::cout << p->get_arrival_time (i) << ",";
        }
#if VERBOSE == 2
        std::cout << " (" << timer.cpu_seconds () << "ms)\n";
#endif
    }

    etavector Vehicle::get_etas ()
    {
        auto stops = _trip->stops ();
        int M (stops.size ());
        etavector etas;
        etas.resize (M);
        if (!valid ()) return etas;
        // std::cout << "\n Vehicle " << _vehicle_id << " =============================";
        std::vector<uint64_t> etam;
        std::vector<double> wts;
        etam.reserve (_state.size ());
        wts.reserve (_state.size ());
        double sumwt;
        for (int i=0; i<M; ++i)
        {
            // need to center each particle's arrival time
            // std::cout << "\n   [" << i << "]: ";
            int tarr = 0;
            int ni = 0;
            etam.clear ();
            wts.clear ();

            for (auto p = _state.begin (); p != _state.end (); ++p)
            {
                // std::cout << p->get_arrival_time (i) << ", ";
                if (p->get_arrival_time (i) > 0)
                {
                    tarr += p->get_weight () * (p->get_arrival_time (i) - _timestamp);
                    etam.push_back (p->get_arrival_time (i));
                    wts.push_back (p->get_weight ());
                }
            }
            if (wts.size () != _N) continue;
            sumwt = std::accumulate (wts.begin (), wts.end (), 0.0);
            if (sumwt < 0.9) continue;

            etas.at (i).stop_id = stops.at (i).stop->stop_id ();
            etas.at (i).estimate = _timestamp + tarr;
            // generate quantiles [0, 5, 50, 95, 100] -> [0, 50, 500, 950, 999]
            std::sort (etam.begin (), etam.end ());
            etas.at (i).quantiles.emplace_back (0.0, etam.front ());
            int qi = 0;
            double cwt = 0.0;
            while (cwt <= 0.05 * sumwt && qi < wts.size ())
            {
                cwt += wts.at (qi);
                qi++;
                // std::cout << ".";
            }
            etas.at (i).quantiles.emplace_back (5.0, etam.at (qi));
            while (cwt <= 0.5 * sumwt && qi < wts.size ())
            {
                cwt += wts.at (qi);
                qi++;
                // std::cout << ".";
            }
            etas.at (i).quantiles.emplace_back (50.0, etam.at (qi));
            while (cwt <= 0.95 * sumwt && qi < wts.size ())
            {
                cwt += wts.at (qi);
                qi++;
                // std::cout << ".";
            }
            etas.at (i).quantiles.emplace_back (95.0, etam.at (qi));
            etas.at (i).quantiles.emplace_back (100.0, etam.back ());
        }
        return etas;
    }

    void Particle::predict_etas (RNG& rng)
    {
        // std::cout << "\n > ";
        if (complete) return;
        // std::cout << "| ";

        // get STOPS
        std::vector<StopTime>* stops;
        if (!vehicle || !vehicle->trip ()) return;
        stops = &(vehicle->trip ()->stops ());
        int M (stops->size ());
        if (M == 0) return;
        at.clear ();
        at.resize (M, 0);
        stop_index = find_stop_index (distance, stops);
        // std::cout << " @" << m << " > ";
        if (stop_index == M-1) 
        {
            return;
        }
        double Dmax = stops->back ().distance;
        
        if (!vehicle->trip ()->shape ()) return;
        std::vector<ShapeSegment>* segments;
        segments = &(vehicle->trip ()->shape ()->segments ());
        int L (segments->size ());
        unsigned int l (find_segment_index (distance, segments));
        double next_segment_d;
        next_segment_d = (l+1 >= L-1) ? Dmax : segments->at (l+1).distance;

        double dcur = distance;
        double dnext;
        uint64_t t0 = vehicle->timestamp ();
        // std::cout << t0 << " > ";
        int etat, tt_total = 0;
        double vel;
        vel = segments->at (l).segment->sample_speed (rng);
        if (vel == 0.0 && speed >= 0.0) vel = speed;
        while (vel <= 0 || vel > 30)
        {
            vel = rng.rnorm () * 8.0 + 15.0;
        }
        while (stop_index < M-1)
        {
            stop_index++; // `next` stop index
            dnext = stops->at (stop_index).distance;
            etat = 0;
            // std::cout << " [";
            while (next_segment_d < dnext && l < L-1)
            {
                // time to get to end of segment
                etat += (next_segment_d - dcur) / vel;
                tt_total += (next_segment_d - dcur) / vel;
                dcur = next_segment_d;
                l++;
                next_segment_d = (l+1 >= L-1) ? Dmax : segments->at (l+1).distance;
                vel = segments->at (l).segment->sample_speed (rng, tt_total);
                if (vel == 0.0 && speed >= 0.0) vel = speed;
                while (vel <= 0.0 || vel > 30)
                {
                    vel = rng.rnorm () * 8.0 + 15.0;
                }
                // std::cout << vel << "; ";
            }
            // std::cout << "] ";
            etat += (dnext - dcur) / vel;
            tt_total += (dnext - dcur) / vel;
            at.at (stop_index) = t0 + etat; // makes no sense because speeds are noise
            dcur = dnext;
            // and add some dwell time
            t0 = at.at (stop_index);
            if (rng.runif () < vehicle->dwell_time ())
            {
                t0 += vehicle->gamma () - vehicle->dwell_time () * log (rng.runif ());
            }
            // std::cout << "(" << m << ") " << at.at (m) << ", ";
        }
    }

}
