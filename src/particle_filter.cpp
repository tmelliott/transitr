/**
 * Define the particle filter functions.
 */
#include "particle_filter.h"
#include "timing.h"

namespace Gtfs {

    void Vehicle::initialize (Event& e, RNG& rng)
    {
        std::cout << "\n    -> initializing";
        initialize (rng);
        _timestamp = e.timestamp;
        _delta = 0;

        // initialize based on the event
        double dmax, dist;
        int x, dr;
        if (e.type == EventType::gps)
        {
            dmax = _trip->shape ()->path ().back ().distance;
            dist = _trip->shape ()->distance_of (_position);
            auto pt = _trip->shape ()->coordinates_of (dist);
            if (distanceEarth (pt, _position) > 50)
            {
                dist = -1.0;
            }
        }
        else
        {
            if (_stop_index >= _trip->stops ().size ()) return;
            dist = _trip->stops ().at (_stop_index).distance;
        }
        
        double d, u;
        Particle* p;
        for (int i=0; i<_N; ++i)
        {
            if (e.type == EventType::gps)
            {
                // initialize each particle within 100m of obs
                u = rng.runif ();
                d = (dist < 0 ? u * dmax : fmin(dmax, u * 200 + dist / 2.0));
                _state.emplace_back (d, rng.runif () * 30.0, rng.rnorm () * _systemnoise, this);
            }
            else
            {
                // stick the bus AT the stop
                _state.emplace_back (dist, rng.runif () * 30, rng.rnorm () * _systemnoise, this);

                // point to the particle
                p = &(_state.back ());
                p->bus_stop (e.timestamp, rng);
                if (e.type == EventType::departure)
                {
                    // shift at <- dt
                    int dwell (p->get_departure_time (e.stop_index) - p->get_arrival_time (e.stop_index));
                    p->set_arrival_time (e.stop_index, p->get_arrival_time (e.stop_index) - dwell);
                    p->set_departure_time (e.stop_index, p->get_departure_time (e.stop_index) - dwell);
                }
                
            }
        }

        double initwt = 1.0 / (double) _N;
        for (auto p = _state.begin (); p != _state.end (); ++p) p->set_weight (initwt);

#if SIMULATION
        {
            double dbar = 0.0, vbar = 0.0;
            for (auto& p : _state)
            {
                dbar += p.get_weight () * p.get_distance ();
                vbar += p.get_weight () * p.get_speed ();
            }

            latlng px = _trip->shape ()->coordinates_of (dbar);

            std::ostringstream fname;
            fname << "history/v" << _vehicle_id << "_mutate.csv";
            std::ofstream fout;
            fout.open (fname.str ().c_str (), std::ofstream::app);
            fout << _timestamp
                << "," << _trip->trip_id ()
                << "," << "initialize"
                << "," << (e.type == EventType::gps ? std::to_string (e.position.latitude) : "")
                << "," << (e.type == EventType::gps ? std::to_string (e.position.longitude) : "")
                << "," << (e.type == EventType::gps ? "" : std::to_string (e.stop_index))
                << "," << dist_to_route
                << "," << _delta
                << "," << dbar
                << "," << vbar
                << "," << (e.type == EventType::gps ? std::to_string (px.latitude) : "")
                << "," << (e.type == EventType::gps ? std::to_string (px.longitude) : "")
                << "," << (e.type == EventType::gps ? std::to_string (distanceEarth (px, _position)) : "")
                << "," << "" // no likelihood
                << "\n";
            fout.close ();
        }
#endif

    }

    void Vehicle::initialize (RNG& rng)
    {
        _state.clear ();
        _state.reserve (_N);
        resample_count = 0;

        if (_trip == nullptr || _trip->shape () == nullptr) return;

        _newtrip = false;
        _complete = false;
        bad_sample = false;
        resample = true;

        _current_segment = 0;
        _current_stop = 0;
        _segment_travel_times.clear ();
        _stop_arrival_times.clear ();

        _segment_travel_times.resize (_trip->shape ()->segments ().size (), 0);
        _stop_arrival_times.resize (_trip->stops ().size (), 0);
    }

    void Vehicle::mutate (RNG& rng, Gtfs* gtfs)
    {
        int nnew = time_events.size () - current_event_index;
        std::cout << "\n\n+ vehicle " << _vehicle_id << ": "       
            << nnew << " new events";
        if (nnew == 0) return;

        std::cout << "\n    [";
        if (_timestamp > 0 && current_event_index > 0) 
        {
            std::cout << _timestamp << "] "
                << _trip->route ()->route_short_name ()
                << " (" << _trip->stops ().at (0).departure_time << "): ";
            time_events.at (current_event_index - 1).print ();
        }
        else std::cout << "uninitialized]";

        // repeat until there are no more events
        while (current_event_index < time_events.size ())
        {
            auto e = time_events.at (current_event_index);

            if (_trip == nullptr || _trip->trip_id () != e.trip_id)
            {
                // assign trip <--> vehicle
                set_trip (gtfs->find_trip (e.trip_id));
                _newtrip = _trip != nullptr;

                _previous_state.clear ();
                _previous_ts = 0;
                _timestamp = 0;
                _state.clear ();
                estimated_dist = 0.0;
            }

            if (_trip == nullptr)
            {
                throw std::runtime_error ("Trip not found");
            }

            std::cout << "\n    [" << e.timestamp << "] "
                << _trip->route ()->route_short_name ()
                << " (" << _trip->stops ().at (0).departure_time << "): ";
            e.print ();

            // is the event "bad"?
            dist_to_route = 0.0;
            switch(e.type)
            {
                case EventType::gps :
                    {
                        // if the distance to the route > 50m, it's bad (skip it)
                        double dist = _trip->shape ()->distance_of (e.position);
                        latlng coord = _trip->shape ()->coordinates_of (dist);
                        dist_to_route = distanceEarth (e.position, coord);
                        std::cout << " -> distance to route = " << dist_to_route << "m";
                        if (dist_to_route > 50)
                        {
                            current_event_index++;
                            return;
                        }
                        break;
                    }
                case EventType::arrival :
                    {
                        // this is tricky ...
                    }
                case EventType::departure :
                    {
                        // no checks
                        break;
                    }
            }


            // pass the event data
            if (e.type == EventType::gps)
            {
                _position = latlng (e.position.latitude, e.position.longitude);
            }
            else
            {
                _stop_index = e.stop_index;
                // ----------------------------------- CHECK STOP INDEX
                std::cout << " (d = " 
                    << _trip->stops ().at (_stop_index).distance << "m)";
            }

            mutate_to (e, rng);

            // if the current iteration fails, start again from here
            if (bad_sample) initialize (e, rng);

            {
                // NETWORK STUFF
                double dmin = _trip->shape ()->path ().back ().distance;
                for (auto p = _state.begin (); p != _state.end (); ++p)
                {
                    if (p->get_distance () < dmin)
                    {
                        dmin = p->get_distance ();
                    }
                }
                std::vector<ShapeSegment>& segs = _trip->shape ()->segments ();
                int m = find_stop_index (dmin, &(_trip->stops ()));
                int l = find_segment_index (dmin, &segs);


                // update segment travel times for intermediate ones ...
                double tt, ttp, err;
                int n;
                while (_current_segment < m)
                {
                    // get the average travel time for particles along that segment
                    tt = 0.0;
                    n = 0;
                    for (auto p = _state.begin (); p != _state.end (); ++p)
                    {
                        ttp = p->get_travel_time (_current_segment) * p->get_weight ();
                        if (ttp > 0)
                        {
                            tt += ttp;
                            n++;
                        }
                    }
                    if (n < _N)
                    {
                        _current_segment++;
                        continue;
                    }

                    err = std::accumulate (_state.begin (), _state.end (), 0.0,
                                           [=](double a, Particle& p) {
                                                return a + p.get_weight () * pow(p.get_travel_time (_current_segment) - tt, 2);
                                           });

                    // if the error is effectively 0 ...
                    if (err < 0.001) err = 10.0;

                    _segment_travel_times.at (_current_segment) = round (tt);
                    segs.at (_current_segment).segment->push_data (tt, err, _timestamp);
                    _current_segment++;
                }
                
                // NOTE: need to ignore segment if previous segment travel time is 0
                // (i.e., can't be sure that the current segment travel time is complete)
                
                
            }

            current_event_index++;
        }

        return;
    }

    void Vehicle::mutate_to (Event& e, RNG& rng)
    {
        if (_newtrip || bad_sample)
        {
            initialize (e, rng);
            return;
        }

        bad_sample = true;

        _delta = e.timestamp - _timestamp;
        _timestamp = e.timestamp;

        // move the particles
        if (_complete || !valid () || _delta == 0) return;
        std::cout << "\n     + " << _delta << " seconds";

        bool all_complete = true;
        double dbar = 0.0, vbar = 0.0, 
            dbar2 = 0.0, vbar2 = 0.0, ddbar = 0.0,
            dtbar = 0;
        for (auto& p : _state)
        {
            if (_N < 20)
                std::cout << "\n      [" << p.get_distance () 
                    << ", " << p.get_speed () 
                    << ", " << (p.get_stop_index () + 1)
                    << "]";
            dbar += p.get_distance () * p.get_weight ();
            vbar += p.get_speed () * p.get_weight ();

            p.travel (_delta, e, rng);

            if (p.is_complete ()) 
            {
                if (_N < 20) std::cout << " -> complete";
                continue;
            }
            // if any aren't complete, prevent vehicle from finishing trip
            all_complete = false;
            if (_N < 20)
                std::cout << " -> [" << p.get_distance () 
                    << ", " << p.get_speed () 
                    << ", " << (p.get_stop_index () + 1)
                    << "]";
            
            dbar2 += p.get_distance () * p.get_weight ();
            vbar2 += p.get_speed () * p.get_weight ();

            // calculate particle likelihood
            switch (e.type)
            {
                case EventType::gps :
                    {
                        double err = fmax (_gpserror, dist_to_route);
                        p.calculate_likelihood (e.position, _trip->shape ()->path (), err);
                        double d = p.get_distance ();
                        auto pos = _trip->shape ()->coordinates_of (d);
                        ddbar += distanceEarth (e.position, pos) * p.get_weight ();
                        break;
                    }
                case EventType::arrival :
                    p.calculate_likelihood (e, _arrival_error);
                    dtbar += (int)(p.get_arrival_time (e.stop_index) - e.timestamp) * p.get_weight ();
                    break;
                case EventType::departure :
                    p.calculate_likelihood (e, _departure_error);
                    dtbar += (int)(p.get_departure_time (e.stop_index) - e.timestamp) * p.get_weight ();
                    break;
            }

            if (_N < 20)
                std::cout << " => l(Y|Xi) = " << exp (p.get_ll ());
        }
        std::cout << "\n    =========================================================================\n"
            << "      [" << dbar << ", " << vbar << "] -> "
            << "[" << dbar2 << ", " << vbar2<< "] => ";
        latlng px = latlng ();
        switch (e.type)
        {
            case EventType::gps :
                {
                    px = _trip->shape ()->coordinates_of (dbar2);
                    std::cout << "d(h(X), y) = " << ddbar 
                        << " [" << px.latitude << ", " << px.longitude << "]";
                    break;
                }
            default :
                std::cout << (e.type == EventType::arrival ? "arrival" : "departure")
                    << " diff " << dtbar << "s";
        }

        _Neff = 0;
        double sumlh = std::accumulate (_state.begin (), _state.end (), 0.0,
                                        [](double a, Particle& p) {
                                            return a + p.get_weight () * exp (p.get_ll ());
                                        });

#if SIMULATION
        {
            std::ostringstream fname;
            fname << "history/v" << _vehicle_id << "_mutate.csv";
            std::ofstream fout;
            fout.open (fname.str ().c_str (), std::ofstream::app);
            fout << _timestamp
                << "," << _trip->trip_id ()
                << "," << e.type_name ()
                << "," << (e.type == EventType::gps ? std::to_string (e.position.latitude) : "")
                << "," << (e.type == EventType::gps ? std::to_string (e.position.longitude) : "")
                << "," << (e.type == EventType::gps ? "" : std::to_string (e.stop_index))
                << "," << dist_to_route
                << "," << _delta
                << "," << dbar2
                << "," << vbar2
                << "," << (e.type == EventType::gps ? std::to_string (px.latitude) : "")
                << "," << (e.type == EventType::gps ? std::to_string (px.longitude) : "")
                << "," << (e.type == EventType::gps ? std::to_string (distanceEarth (px, _position)) : "")
                << "," << sumlh
                << "\n";
            fout.close ();
        }
#endif

        std::cout << "\n   -> sum(l(y|x)) = " << sumlh;
        // if no likelihoods are that big, give up
        if (sumlh < 1e-16)
        {
            if (n_bad < 2)
            {
                // we effectively want to ignore that observation
                bad_sample = false;
                n_bad++;   
            }
            return;
        }
        else
        {
            n_bad = 0;
        }

        // normalize weights
        std::vector<double> wt;
        wt.reserve (_state.size ());
        for (auto& p : _state)
        {
            p.set_weight (p.get_weight () * exp (p.get_ll ()) / sumlh);
            wt.push_back (p.get_weight ());
        }
        double sumwt = std::accumulate (wt.begin (), wt.end (), 0.0);

        // print posterior values
        dbar = 0.0;
        vbar = 0.0;
        ddbar = 0.0;
        dtbar = 0;
        for (auto& p : _state)
        {
            dbar += p.get_distance () * p.get_weight ();
            vbar += p.get_speed () * p.get_weight ();
            switch (e.type)
            {
                case EventType::gps :
                    {
                        double d = p.get_distance ();
                        auto pos = _trip->shape ()->coordinates_of (d);
                        ddbar += distanceEarth (e.position, pos) * p.get_weight ();
                        break;
                    }
                case EventType::arrival :
                    dtbar += (int)(p.get_arrival_time (e.stop_index) - e.timestamp) * p.get_weight ();
                    break;
                case EventType::departure :
                    dtbar += (int)(p.get_departure_time (e.stop_index) - e.timestamp) * p.get_weight ();
                    break;
            }
        }
        std::cout << "\n         => Posterior: ["
            << dbar << ", " << vbar << "] => ";
        switch (e.type)
        {
            case EventType::gps :
                {
                    std::cout << "(getpx_len="
                        << _trip->shape ()->path ().size () << ")";
                    px = _trip->shape ()->coordinates_of (dbar);
                    std::cout << "got";
                    std::cout.flush ();
                    std::cout << "d(h(X), y) = " << ddbar 
                        << " [" << px.latitude << ", " << px.longitude << "]";
                    break;
                }
            default :
                std::cout << (e.type == EventType::arrival ? "arrival" : "departure")
                    << " diff " << dtbar << "s";
        }

        // sum(wt) in (0.999, 1.0001)
        std::cout << "\n   -> sum(wt) = " << sumwt;
        if (fabs (1 - sumwt) > 1e-4) return;
        bad_sample = false;

        // calcualte Neff.
        double sumwt2 = std::accumulate(wt.begin (), wt.end (), 0.0,
                                        [](double a, double b) {
                                            return a + pow(b, 2);
                                        });

        _Neff = pow (sumwt2, -1);
        std::cout << "\n   -> Neff = " << _Neff;

#if SIMULATION
        {
            std::ostringstream fname;
            fname << "history/v" << _vehicle_id << "_update.csv";
            std::ofstream fout;
            fout.open (fname.str ().c_str (), std::ofstream::app);
            fout << _timestamp
                << "," << _trip->trip_id ()
                << "," << e.type_name ()
                << "," << dbar
                << "," << vbar
                << "," << (e.type == EventType::gps ? std::to_string (px.latitude) : "")
                << "," << (e.type == EventType::gps ? std::to_string (px.longitude) : "")
                << "," << (e.type == EventType::gps ? std::to_string (distanceEarth (px, _position)) : "")
                << "," << _Neff
                << "," << (_Neff >= (_N / 4) ? 0 : 1)
                << "\n";
            fout.close ();
        }
#endif


        if (_Neff >= (_N / 4)) return;

        std::cout << " -> resampling";
        select (rng);

    }

    void Vehicle::select (RNG& rng)
    {
        std::vector<double> wt;
        wt.reserve (_state.size ());
        double sumwt = 0.0;
        for (auto p = _state.begin (); p != _state.end (); ++p)
        {
            wt.push_back (p->get_weight ());
            sumwt += p->get_weight ();
        }

        // resample u \in [0, 1)
        double u;
        unsigned int j;
        std::vector<Particle> _newstate;
        _newstate.reserve (_N);
        for (int i=0; i<_N; ++i)
        {
            u = rng.runif () * sumwt;
            for (j=0; j<_N; ++j)
            {
                if (u < wt.at (j))
                {
                    _newstate.emplace_back (_state.at (j));
                    break;
                }
                // if not less, subtract wt from u
                u -= wt.at (j);
            }
        }

        _state.clear ();
        _state = std::move(_newstate);

        double initwt = pow(_N, -1);
        for (auto& p : _state) p.set_weight (initwt);
        resample_count++;
    }

    void Vehicle::predict_etas (RNG& rng)
    {
        if (!valid ()) return;

#if VERBOSE == 2
        Timer timer;
        std::cout << "- vehicle " << _vehicle_id << " - predicting etas";
#endif
        for (auto p = _state.begin (); p != _state.end (); ++p) 
        {
            p->predict_etas (rng);
            // std::cout << "\n";
            // for (auto eta : p->get_arrival_times ()) std::cout << eta << ", ";
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

    double Vehicle::distance ()
    {
        double distance = 0.0;
        for (auto p = _state.begin (); p != _state.end (); ++p) 
            distance += p->get_weight () * p->get_distance ();
        return distance;
    }

    double Vehicle::speed ()
    {
        double speed = 0.0;
        for (auto p = _state.begin (); p != _state.end (); ++p) 
            speed += p->get_weight () * p->get_speed ();
        return speed;
    }

    int Vehicle::progress ()
    {
        double d = distance ();
        double dmax = _trip->shape ()->path ().back ().distance;
        return (100 * d / dmax) + 0.5;
    }



    void Particle::travel (int delta, Event& e, RNG& rng)
    {
        if (complete || !vehicle || !vehicle->trip () || 
            !vehicle->trip ()->shape ()) return;

        // do the particle physics
        double Dmax = vehicle->trip ()->shape ()->path ().back ().distance;
        if (distance >= Dmax) 
        {
            distance = Dmax;
            complete = true;
            return;
        }
        
        // get STOPS
        std::vector<StopTime>* stops;
        stops = &(vehicle->trip ()->stops ());
        int M (stops->size ());
        // if (stop_index == 0) stop_index = find_stop_index (distance, stops);
        // std::cout << " [M=" << M << ",m=" << m << "], ";
        if (stop_index == M-1) 
        {
            distance = Dmax;
            complete = true;
            return;
        }

        double next_stop_d = stops->at (stop_index + 1).distance;
        
        // get SEGMENTS
        std::vector<ShapeSegment>* segments;
        segments = &(vehicle->trip ()->shape ()->segments ());
        int L (segments->size ());
        unsigned int l (find_segment_index (distance, segments));
        // std::cout << " [L=" << L << ",l=" << l << "], ";
        double next_segment_d;
        next_segment_d = (l+1 >= L-1) ? Dmax : segments->at (l+1).distance;

        // std::cout << " tt.size = " << tt.size ()
            // << ", at.size = " << at.size ();
        
        // allow vehicle to remain stationary if at a stop:
        // if (distance == stops->at (stop_index).distance)
        // {
        //     if (rng.runif () < 0.05)
        //     {
        //         double w = vehicle->gamma () - vehicle->dwell_time () * log (rng.runif ());
        //         delta = fmax (0, delta - round (w));
        //         // we don't want this to affect the speed
        //     }
        // }
        // else if (distance == segments->at (l).distance &&
        //          rng.runif () < 0.05)
        // {
        //     double w = - log (rng.runif ()) * delta;
        //     delta = fmax (0, delta - round (w));
        //     if (tt.at (l) >= 0)
        //         tt.at (l) = tt.at (l) + 1;
        // }

         
        // adjust "delta" if arrival/departure already set
        uint64_t tx = vehicle->timestamp () - delta;
        if (dt.at (stop_index) > tx)
        {
            // particle is AT the stop, not yet departed (it's just about to...)
            if (e.type == EventType::departure && stop_index == e.stop_index)
            {
                uint64_t newtime = rng.rnorm () * vehicle->departure_error () + e.timestamp + 0.5;
                // departure time has to be >= arrival time
                dt.at (stop_index) = fmax (newtime, at.at (stop_index));
                delta = (int)(vehicle->timestamp () - dt.at (stop_index));
            }
            // can't fudge arrival time (in case its based on speed!!!)
        }
        else if (at.at (stop_index) > tx)
        {
            // particle has not yet ARRIVED at the stop
            if (e.type == EventType::departure && stop_index == e.stop_index)
            {
                // let's fudge 
                uint64_t newtime = rng.rnorm () * vehicle->departure_error () + e.timestamp + 0.5;
                dt.at (stop_index) = fmax (newtime, at.at (stop_index));
                delta = (int)(vehicle->timestamp () - dt.at (stop_index));
            }
            else if (e.type == EventType::arrival && stop_index == e.stop_index)
            {
                // resample dwell time (to allow variation)
                bus_stop (at.at (stop_index), rng);
                delta = (int)(vehicle->timestamp () - dt.at (stop_index));
            }
            else
            {
                delta = (int)(vehicle->timestamp () - at.at (stop_index));
            }
        }

        if (distance > stops->at (stop_index).distance && rng.runif () < 0.01)
        {
            // a very small chance for particles to remain stationary
            int tr = rng.runif () * delta + 0.5;
            if (tt.at (l) >= 0)
                tt.at (l) = tt.at (l) + tr;
            delta -= tr;

            // speed = 0.0;
            // when /not/ at a bus stop, set speed to 0 and wait
            // double w = - log (rng.runif ()) * delta;
            // delta = fmax (0.0, delta - round (w));
            // then the bus needs to accelerate back up to speed ... for how many seconds?
            // accelerating = 5.0 + rng.runif () * 10.0;
            // acceleration = 2.0 + rng.rnorm () * vehicle->system_noise ();
        }


        double speed_mean = 10.0;
        double speed_sd = 100.0;
        // if (segments->at (l).segment->travel_time () > 0 &&
        //     segments->at (l).segment->uncertainty () > 0) 
        // {
        //     speed_mean = segments->at (l).segment->get_speed ();
        //     speed_sd = - segments->at (l).segment->length () / pow (speed_mean, 2) *
        //         segments->at (l).segment->uncertainty ();
        //     speed_sd = pow(speed_sd, 0.5);
        // }
        double vmax = 30; //rng.runif () < 0.5 ? 30.0 : 15.0;

        // while (distance < Dmax && delta > 0.0)
        while (behind_event (e, delta))
        {
            // add system noise to acceleration to ensure speed remains in [0, vmax]
            double accel_prop (-100.0);
            double n = 0;
            if (speed > vmax)
            {
                speed = rng.runif () * vmax;
            }
            while (speed + accel_prop < 0.0 || speed + accel_prop > vmax && n < 1000)
            {
                accel_prop = rng.rnorm () * vehicle->system_noise () * 
                    (1.0 + (double)n / 100.0);
                n++;
                // if (accelerating > 0.0)
                // {
                //     accel_prop += acceleration;
                //     accelerating--;
                // }
            }

            // double v = fmax (0, fmin (30, speed + acceleration));
            // double v = speed;
            // double vstar = speed + accel_prop;
            // double alpha = (pow (v - speed_mean, 2) - pow(vstar - speed_mean, 2)) / (2 * pow (speed_sd, 2));
            // alpha = fmin (0, alpha);
            // if (rng.runif () < exp (alpha))
            // {
            //     acceleration = accel_prop;
            //     speed = vstar;
            // }
            // else
            // {
            //     speed = v;
            // }

            speed += accel_prop;
            distance += speed;
            if (distance >= Dmax) complete = true;
            delta--;

            if (tt.at (l) >= 0)
                tt.at (l) = tt.at (l) + 1;
            
            if (l < L-1 && distance >= next_segment_d)
            {
                // reaching intersection ... 
                l++;
                tt.at (l) = 0;
                next_segment_d = (l+1 >= L-1) ? Dmax : segments->at (l+1).distance;
                // vmax = rng.runif () < 0.5 ? 30.0 : 15.0;
                // if (segments->at (l).segment->travel_time () > 0 &&
                //     segments->at (l).segment->uncertainty () > 0) 
                // {
                //     speed_mean = segments->at (l).segment->get_speed ();
                //     speed_sd = - segments->at (l).segment->length () / pow (speed_mean, 2) *
                //         segments->at (l).segment->uncertainty ();
                //     speed_sd = pow(speed_sd, 0.5);
                // }
                // else
                // {
                    speed_mean = 10.0;
                    speed_sd = 100.0;
                // }
            }

            if (distance >= next_stop_d)
            {
                stop_index++;
                // bus arrives at stop - either stops at stop, of drives past
                if (bus_stop (vehicle->timestamp () - delta, rng)) 
                    distance = next_stop_d;
                // end of route? end.
                if (stop_index >= M-1) break;
                // else update next stop and delta
                next_stop_d = stops->at (stop_index + 1).distance;
                delta -= (int)(dt.at (stop_index) - at.at (stop_index));
                continue;
            }
        }

        // almost at next stop ...
        if (stop_index < M-1 &&
            next_stop_d - distance < 20)
        {
            stop_index++;
            bus_stop (vehicle->timestamp () - delta, rng);
            distance = next_stop_d;
            if (stop_index < M - 1) 
                next_stop_d = stops->at (stop_index + 1).distance;
        }
    }

    /**
     * Simulate the behaviour of a bus when it reaches a bus stop (at time t)
     * @param stop_index the 0-based index of the stop
     * @param time       the UNIX time the particle arrives
     * @param rng        RNG
     */
    bool Particle::bus_stop (uint64_t time, RNG& rng)
    {
        if (stop_index >= at.size ()) return true;

        // bus is arriving
        at.at (stop_index) = time;

        double u (rng.runif ());
        if (u < vehicle->pr_stop ())
        {
            // bus stops
            double dwell = vehicle->gamma () - vehicle->dwell_time () * log (rng.runif ());
            dt.at (stop_index) = time + round (dwell);
            return true;
        }
        else
        {
            // bus doesn't stop
            dt.at (stop_index) = time;
            return false;
        }
    }

    bool Particle::behind_event (Event& e, double delta)
    {
        if (delta > 0) return true;
        if (complete)
        {
            return false;
        }

        switch (e.type)
        {
            case EventType::arrival :
                if (e.stop_index >= at.size ()) return false;
                if (e.stop_index < stop_index) return false;
                return at.at (e.stop_index) == 0;

            case EventType::departure :
                if (e.stop_index >= dt.size ()) return false;
                if (e.stop_index < stop_index) return false;
                return dt.at (e.stop_index) == 0;

            default:
                return false;
        }
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

    void
    Particle::calculate_likelihood (latlng& y, 
                                    std::vector<ShapePt>& path, 
                                    double sigma)
    {
        latlng ppos = vehicle->trip ()->shape ()->coordinates_of (distance);
        // (log) distance between points
        double ld = log (distanceEarth (ppos, vehicle->position ()));
        // (log) (d/sigma)^2 ~ Chi2(2) ~ Exp(2)
        double lX2 = 2 * (ld - log (sigma));
        // log pdf of lX2 ~ Exp(2)
        log_likelihood = log (0.5) - 0.5 * exp(lX2);

        if (vehicle->get_n () < 20) std::cout << " => d(h(x), y) = " << exp (ld) << "m";
    }

    void Particle::calculate_likelihood (Event& e, double error)
    {
        uint64_t t = (e.type == EventType::arrival) ? 
            get_arrival_time (e.stop_index) :
            get_departure_time (e.stop_index);

        if (t == 0) 
        {
            if (vehicle->get_n () < 20)
                std::cout << " => hasn't "
                    << (e.type == EventType::arrival ? "arrived" : "departed");
            // log_likelihood = 0.0;
            return;
        }

        int tdiff = t - e.timestamp;
        if (vehicle->get_n () < 20)
            std::cout << " => "
                << (e.type == EventType::arrival ? "arrived" : "departed")
                << " at " << t
                << " (diff: " << (tdiff) << "s)";

        log_likelihood = - 0.5 * log (2 * M_PI) - log (error) - 
            pow (tdiff, 2) / 2 / pow (error, 2);
        
    }

    void Particle::calculate_arrival_likelihood (int index, uint64_t time, double error)
    {
        // particle's arrival time at stop INDEX
        if (get_arrival_time (index) == 0) log_likelihood = 0;
        log_likelihood = - 0.5 * log (2 * M_PI) - log (error) - pow (get_arrival_time (index) - time, 2) / 2 / pow (error, 2);
    }
    void Particle::calculate_departure_likelihood (int index, uint64_t time, double error)
    {
        if (get_departure_time (index) == 0) log_likelihood = 0;
        log_likelihood = - 0.5 * log (2 * M_PI) - log (error) - pow (get_departure_time (index) - time, 2) / 2 / pow (error, 2);
    }

    void Particle::set_weight (double w)
    {
        weight = w;
    }

}; // namespace Gtfs
