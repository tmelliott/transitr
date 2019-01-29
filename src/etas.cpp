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
        if (!valid () || complete ()) return;


        Timer timer;
        std::cout << "\n- vehicle " << _vehicle_id << " - predicting etas";

        // std::cout << "\n\n Segment information ...";
        auto segments = _trip->shape ()->segments ();
        // for (int l = 0; l < segments.size (); l++)
        // {
        //     std::cout << "\n  - " << l << ":"
        //         // << " distance = " << segments.at (l).distance
        //         // << ", length = " << segments.at (l).segment->length () << ", "
        //         << " tt=" << segments.at (l).segment->travel_time ()
        //         << " (" << segments.at (l).segment->uncertainty ()
        //         << ")";
        // }

        auto stops = _trip->stops ();
        // std::cout << "\n\n Stop information ...";
        // for (auto j = 0; j < stops.size (); j++)
        // {
        //     std::cout << "\n  - " << j
        //         << ": distance = " << stops.at (j).distance;
        // }

        std::cout << "\n    --- Particle ETAs ...";
        for (auto p = _state.begin (); p != _state.end (); ++p) 
        {
            p->predict_etas (rng);
            // std::cout << "\n  --- stops:";
            // for (int l = 0; l < stops.size (); l++)
            // {
            //     std::cout << "\n     [" << l << "]:"
            //         << p->get_arrival_time (l) << ", "
            //         << p->get_departure_time (l);
            // }

            std::cout << "\n    => ";
            for (auto eta : p->get_arrival_times ()) {
                if (eta <= _timestamp) std::cout << "*, ";
                else std::cout << ((eta - _timestamp)) << ", ";
            }
        }

        std::cout << "\n\n    --- Particle travel times ...";
        // get travel time for each stop
        int L = _trip->shape ()->segments ().size ();
        for (auto& p : _state)
        {
            std::cout << "\n    => ";
            for (int l=0; l<L; l++)
            {
                std::cout << p.get_travel_time_prediction (l) << ", ";
            }
        }

        /**
         * Now, we assume the particles have taken into account any correlation
         * structure between segment travel times.
         *
         * Simply need to compute B and Cov matrix for travel times (for this vehicle)
         */
        
        // auto stops = _trip->stops ();
        int M = stops.size ();
        
        std::cout << std::setprecision (0) << std::fixed;
        std::cout << "\n\n  >> E(B) vector: [";
        for (auto b : _tt_state) std::cout << " " << b << " ";
        std::cout << "]\n  >> Var(B) matrix: ";
        for (auto br : _tt_cov)
        {
            std::cout << "\n  ";
            for (auto bc : br) std::cout << std::setw(9) << std::round (bc) << "  ";
        }
        
        std::cout << "\n\n  >> generate observation of B (stop-stop travel times) and estimate of R\n";
        std::vector<double> tt_obs (M, 0.0);
        std::vector<std::vector<double> > tt_r (M, std::vector<double> (M, 0.0));

        double cov_lj;
        std::cout << "  > curtime = " << _timestamp << "\n";
        std::cout << "  > curstop = " << _current_stop << "\n";
        for (int l=0; l<M; l++)
        {
            std::cout << "\n STOP " << l << ": ";
            tt_obs.at (l) = std::accumulate(_state.begin (), _state.end (), 0.0,
                                            [&](double d, Particle& p) {
                                                int dt;
                                                // std::cout << "[" << p.get_stop_index () << "]";
                                                if (l <= p.get_stop_index ())
                                                {
                                                    // this stop is behind this particle
                                                    dt = 0;
                                                }
                                                if (l == p.get_stop_index () + 1)
                                                {
                                                    // it's the next stop for this particle
                                                    dt = p.get_arrival_time (l) - _timestamp;
                                                }
                                                else 
                                                {
                                                    dt = p.get_arrival_time (l) - p.get_departure_time (l - 1);
                                                    // std::cout << "("
                                                    //     << p.get_arrival_time (l) << " - "
                                                    //     << p.get_departure_time (l - 1) << " = ) ";
                                                }
                                                if (dt < 0) dt = 0;
                                                std::cout << std::setw (5) << dt << ", ";
                                                return d + dt;
                                            });
            tt_obs.at (l) /= (double)_state.size ();
        }
        std::cout << "\n\n Z vector: [";
        for (auto z : tt_obs) std::cout << " " << z << " ";

        for (int l=0; l<M; l++) {
            for (int j=0; j<M; j++)
            {
                cov_lj = std::accumulate(_state.begin (), _state.end (), 0.0,
                                         [&](double d, Particle& p) {
                                            if (l <= p.get_stop_index () || 
                                                j <= p.get_stop_index () ||
                                                p.get_arrival_time (l) == 0) return d;
                                            
                                            double xi, yi;
                                            if (l == p.get_stop_index () + 1)
                                            {
                                                xi = p.get_arrival_time (l) - _timestamp;
                                            }
                                            else
                                            {
                                                xi = p.get_arrival_time (l) - p.get_departure_time (l - 1);
                                            }
                                            
                                            if (l == j) 
                                            {
                                                yi = xi;
                                            }
                                            else if (j == p.get_stop_index () + 1)
                                            {
                                                yi = p.get_arrival_time (j) - _timestamp;
                                            }
                                            else
                                            {
                                                yi = p.get_arrival_time (j) - p.get_departure_time (j - 1);
                                            }

                                            xi -= tt_obs.at (l); // xi - xbar
                                            yi -= tt_obs.at (j); // yi - ybar
                                            return d + xi * yi;
                                         });
                tt_r.at (l).at (j) = cov_lj / (_state.size () - 1);
            }
        }

        std::cout << "]\n R matrix: ";
        for (auto br : tt_r)
        {
            std::cout << "\n  ";
            for (auto bc : br) std::cout << std::setw(9) << std::round (bc) << "  ";
        }
        std::cout << "\n ---------\n";
        

        std::cout << "\n Now update the travel time estimates' state ...\n";
        if (_tt_time == 0)
        {
            _tt_state = tt_obs;
            _tt_cov = tt_r;
        }
        else
        {
            int tt_delta = _timestamp - _tt_time;
            std::cout << " -> " << tt_delta << " seconds ... \n";
            // construct F matrix
            std::cout << "\n  >> F (transition) matrix:";
            
        }


        _tt_time = _timestamp;
        std::cout << std::setprecision (6);

#if VERBOSE == 2
        std::cout << "\n   (" << timer.cpu_seconds () << "ms)\n";
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
        
        // predict travel times along all remaining segments
        if (!vehicle->trip ()->shape ()) return;
        std::vector<ShapeSegment>* segments;
        segments = &(vehicle->trip ()->shape ()->segments ());
        int L (segments->size ());
        ttpred.resize (L, 0);

        // current segment is partial
        // std::cout << "\n\n   --- step 1: travel times from segment " << segment_index << "\n";
        ttpred.at (segment_index) = (distance - segments->at (segment_index).distance) / speed + 0.5;
        // store cumlative travel time for forecasting ahead
        int tcum = ttpred.at (segment_index);
        for (int l=segment_index+1; l<L; l++)
        {
            ttpred.at (l) = segments->at (l).segment->sample_travel_time (rng, tcum);
            if (ttpred.at (l) == 0)
            {
                // something from [5, 25]
                ttpred.at (l) = segments->at (l).segment->length () / (rng.runif () * 20.0 + 15.0);
            }
            // std::cout << "   [" << std::setw(2) << l << "] "
            //     << std::setw (5) << ttpred.at (l);
            tcum += ttpred.at (l);
        }

        // now use those predictions for stop arrival times
        std::vector<StopTime>* stops;
        if (!vehicle || !vehicle->trip ()) return;
        stops = &(vehicle->trip ()->stops ());
        int M (stops->size ());
        if (M == 0) return;
        at.resize (M, 0);
        dt.resize (M, 0);

        double dcur = distance;

        int si (stop_index + 1); // next stop to predict
        unsigned int l, li; // the segment index of current/next stop, respectively        
        double Dmax = stops->back ().distance;

        li = find_segment_index (stops->at (stop_index).distance, segments);

        uint64_t t0 = vehicle->timestamp ();
        int etat; // the eta (in seconds)
        // if (vehicle_at_stop) then add dwell time for stop (stop_index)
        
        double v;
        // std::cout << "\n\n   --- step 2: stop arrival(/departure) predictions (last stop = " 
        //     << stop_index << ")\n";
        // std::cout << " (start at " << t0 << ")";
        while (si < M)
        {
            etat = 0;
            l = li;
            li = find_segment_index (stops->at (si).distance, segments);

            // std::cout << "\n   [" << std::setw (2) << si << "]";
            if (l == li)
            {
                // std::cout << "A ";
                // both stops in the same segment
                v = (l == segment_index) ? speed : segments->at (l).segment->length () / ttpred.at (l);
                etat = (stops->at (si).distance - dcur) / v + 0.5;
                // std::cout << (stops->at (si).distance - dcur) << "m @ "
                //     << v << "m/s -> "
                //     << etat << "s";
            }
            else 
            {
                // std::cout << "B ["
                //     << l << "-" << li << "] ";
                // stops in different segments
                // first, rest of segment (l)
                v = (l == segment_index) ? speed : segments->at (l).segment->length () / ttpred.at (l);
                // using segment length avoids l+1 being out of index
                etat += (segments->at (l).distance + segments->at (l).segment->length () - dcur) / v + 0.5;

                // std::cout << (segments->at (l).distance + segments->at (l).segment->length () - dcur)
                //     << "m @ " << v << "m/s -> "
                //     << etat << "s";
                // then any intermediate segments
                for (int j=l+1; j<li; j++)
                {
                    etat += ttpred.at (j);
                    // std::cout << " + " << ttpred.at (j);
                }

                // then beginning bit of segment (li)
                // if (segments->at (li).distance < stops->at (si).distance)
                // {
                // }
                v = segments->at (li).segment->length () / ttpred.at (li);
                etat += (stops->at (si).distance - segments->at (li).distance) / v + 0.5;
                // std::cout << " + "
                //     << (stops->at (si).distance - segments->at (li).distance)
                //     << "m @ " << v << "m/s -> "
                //     << (stops->at (si).distance - segments->at (li).distance) / v + 0.5
                //     << "s";
            }

            dcur = stops->at (si).distance;
            t0 += etat;
            at.at (si) = t0;
            // std::cout << " = " << etat << " -> " << t0;
            // add dwell time 
            t0 += 0;
            dt.at (si) = t0;
            si++;
        }

        return;

        
        // double dcur = distance;
        // double dnext;
        // uint64_t t0 = vehicle->timestamp ();
        // // std::cout << t0 << " > ";
        // int etat, tt_total = 0;
        // // double vel;
        // // vel = segments->at (l).segment->sample_speed (rng);
        // // if (vel == 0.0 && speed >= 0.0) vel = speed;
        // // // std::cout << "\n        P: stop " << (stop_index+1) << "/" << M
        // // //     << ", segment " << (l+1) << "/" << L;
        // // while (vel <= 0 || vel > 30)
        // // {
        // //     vel = rng.rnorm () * 8.0 + 15.0;
        // // }

        // // -------------------------------- FOR NOW! update in future
        // ttpred.resize (L, 0);

        // std::cout << "\n";
        // while (si < M-1)
        // {
        //     std::cout << " | si = " << si << " -> ";
        //     si++; // `next` stop index
        //     std::cout << si << ", ";
        //     dnext = stops->at (si).distance;
        //     etat = 0;
        //     // std::cout << " [";
        //     while (next_segment_d <= dnext && l < L-1)
        //     {
        //         std::cout << " l = " << l << " -> ";
        //         // time to get to end of segment
        //         ttpred.at (l) += (next_segment_d - dcur) / vel + 0.5; // 0.5 for integer rounding
        //         etat += ttpred.at (l);
        //         tt_total += ttpred.at (l);
        //         dcur = next_segment_d;
        //         l++;
        //         std::cout << l << ", ";
        //         next_segment_d = (l+1 >= L-1) ? Dmax : segments->at (l+1).distance;
        //         vel = segments->at (l).segment->sample_speed (rng, tt_total);
        //         if (vel == 0.0 && speed >= 0.0) vel = speed;
        //         while (vel <= 0.0 || vel > 30)
        //         {
        //             vel = rng.rnorm () * 8.0 + 15.0;
        //         }
        //         // std::cout << vel << "; ";
        //     }
        //     // std::cout << "] ";
        //     ttpred.at (l) += (dnext - dcur) / vel + 0.5;
        //     etat += (dnext - dcur) / vel + 0.5;
        //     tt_total += (dnext - dcur) / vel + 0.5;
        //     t0 += etat;
        //     dcur = dnext;
        //     at.at (si) = t0;
        //     // and add some dwell time
        //     if (rng.runif () < vehicle->pr_stop ())
        //     {
        //         t0 += vehicle->gamma () - vehicle->dwell_time () * log (rng.runif ());
        //     }
        //     dt.at (si) = t0;
        //     // std::cout << "(" << m << ") " << at.at (m) << ", ";
        // }
        // // last segment
        // std::cout << l << "..." << ttpred.at (l) << " ... ";
        // if (l < L-1)
        // {
        //     std::cout << " not finished...";
        // }
    }

}
