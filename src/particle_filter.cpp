/**
 * Define the particle filter functions.
 */
#include "particle_filter.h"
#include "timing.h"

namespace Gtfs {

    void Vehicle::initialize (Event& e, RNG& rng)
    {
        initialize (rng);
        _timestamp = e.timestamp;

        // initialize based on the event
        double dmax, dist;
        if (e.type == EventType::gps)
        {
            dmax = _trip->shape ()->path ().back ().distance;
            dist = _trip->shape ()->distance_of (_position);
            auto pt = _trip->shape ()->coordinates_of (dist);
            if (distanceEarth (pt, _position) > 50)
            {
                dist = 0.0;
            }
        }
        else
        {
            if (_stop_index >= _trip->stops ().size ()) return;
            dist = _trip->stops ().at (_stop_index).distance;
        }
        
        for (int i=0; i<_N; ++i)
        {
            switch (e.type)
            {
                case EventType::gps :
                    // initialize each particle within 100m of obs
                    double d, u;
                    u = rng.runif ();
                    d = (dist == 0 ? u * dmax : fmin(dmax, u * 200 + dist / 2.0));
                    _state.emplace_back (d, rng.runif () * 30.0, rng.rnorm () * _systemnoise, this);
                    break;
                case EventType::arrival :
                case EventType::departure :
                    // initialize points at the stop, I guess ...
                    _state.emplace_back (dist, rng.runif () * 30, rng.rnorm () * _systemnoise, this);
                    break;
            }
        }

        double initwt = 1.0 / (double) _N;
        for (auto p = _state.begin (); p != _state.end (); ++p) p->set_weight (initwt);

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
        if (_timestamp > 0) 
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

            // pass the event data
            if (e.type == EventType::gps)
            {
                _position = latlng (e.position.latitude, e.position.longitude);
            }
            else
            {
                _stop_index = e.stop_index;
            }

            std::cout << "\n    [" << e.timestamp << "] "
                << _trip->route ()->route_short_name ()
                << " (" << _trip->stops ().at (0).departure_time << "): ";
            e.print ();

            if (e.type != EventType::gps)
            {
                std::cout << " (d = " 
                    << _trip->stops ().at (_stop_index).distance << "m)";
            }

            mutate_to (e, rng);

            current_event_index++;
        }

        return;
        // // are there any stop updates that need accounting for?
        // if (_stop_time_updates.size () > 0 && 
        //     _last_stop_update_index >= 0 && 
        //     _last_stop_update_index < _stop_time_updates.size ())
        // {
        //     STU& stu = _stop_time_updates.at (_last_stop_update_index);
        //     if (stu.timestamp > 0)
        //     {
        //         // most recent update (and maybe others!!) needs to be incorporated into the likelihood
        //         // BUT--- for now, just using a DIFFERENT likelihood for these
        //         if (!stu.used_arrival && stu.arrival_time > 0)
        //         {
        //             _skip_observation = true;
        //             if (_timestamp == stu.arrival_time)
        //             {
        //                 // Y is from arrival, so forget it  - timestamp and delta are correct
        //                 std::cout << " - Arr1 d=" << _delta << "\n";
        //                 mutate2 (rng);
        //                 _delta = 0;
        //             }
        //             else if (_timestamp < stu.arrival_time)
        //             {
        //                 std::cout << " - Arr2 d=" << _delta << "\n";
        //                 // mutate to Y, then arrival
        //                 mutate2 (rng);
        //                 _delta = stu.arrival_time - _timestamp;
        //                 _previous_ts = _timestamp;
        //                 _timestamp = stu.arrival_time;
        //                 mutate2 (rng);
        //                 _delta = 0;
        //             }
        //             else
        //             {
        //                 // mutate to arrival, then Y
        //                 auto ts = _timestamp;
        //                 _delta = stu.arrival_time - _previous_ts;
        //                 _timestamp = stu.arrival_time;
        //                 std::cout << " - Arr3 d=" << _delta << "\n";
        //                 mutate2 (rng);
        //                 _timestamp = ts;
        //                 _delta = _timestamp - stu.arrival_time;
        //                 _previous_ts = stu.arrival_time;
        //                 mutate2 (rng);
        //                 _delta = 0;
        //             }
        //             stu.used_arrival = true;
        //         }

        //         if (!stu.used_departure && stu.departure_time > 0)
        //         {
        //             _skip_observation = true;
        //             if (_timestamp == stu.departure_time)
        //             {
        //                 std::cout << " - Dep1 d=" << _delta << "\n";
        //                 // Y is from departure, so forget it - timestamp and delta are correct
        //                 mutate2 (rng);
        //             }
        //             else if (_timestamp < stu.departure_time)
        //             {
        //                 // mutate to Y, then departure
        //                 std::cout << " - Dep2 d=" << _delta << "\n";
        //                 if (_delta > 0)
        //                 {
        //                     mutate2 (rng);
        //                 }
        //                 _delta = stu.departure_time - _timestamp;
        //                 _previous_ts = _timestamp;
        //                 _timestamp = stu.departure_time;
        //                 mutate2 (rng);
        //                 _delta = 0;
        //             }
        //             else
        //             {
        //                 auto ts = _timestamp;
        //                 _delta = stu.departure_time - _previous_ts;
        //                 _timestamp = stu.departure_time;
        //                 std::cout << " - Dep3 d=" << _delta << "\n";
        //                 mutate2 (rng);
        //                 _timestamp = ts;
        //                 _delta = _timestamp - stu.departure_time;
        //                 _previous_ts = stu.departure_time;
        //                 mutate2 (rng);
        //                 _delta = 0;
        //             }
        //             stu.used_departure = true;
        //         }
        //     }
        // }

        // // ok, there was no arrival stuff 
        // if (!_skip_observation) 
        // {
        //     std::cout << " - Obs d=" << _delta << "\n";
        //     mutate2 (rng);
        // }

        // _skip_observation = false;
    }

    void Vehicle::mutate_to (Event& e, RNG& rng)
    {
        if (_newtrip || bad_sample)
        {
            _delta = 0;
            std::cout << "\n    -> initializing";
            initialize (e, rng);
            return;
        }

        _delta = e.timestamp - _timestamp;
        _timestamp = e.timestamp;

        // move the particles
        if (_complete || !valid () || _delta == 0) return;
        std::cout << "\n     + " << _delta << " seconds";

        bool all_complete = true;
        for (auto& p : _state)
        {
            std::cout << "\n      [" << p.get_distance () 
                << ", " << p.get_speed () 
                << ", " << p.get_stop_index ()
                << "]";
            p.travel (_delta, rng);
            if (p.is_complete ()) 
            {
                std::cout << " -> complete";
                continue;
            }
            // if any aren't complete, prevent vehicle from finishing trip
            all_complete = false;
            std::cout << " -> [" << p.get_distance () 
                << ", " << p.get_speed () 
                << ", " << p.get_stop_index ()
                << "]";

            // calculate particle likelihood
            switch (e.type)
            {
                case EventType::gps :
                    p.calculate_likelihood (e.position, _trip->shape ()->path (), _gpserror);
                    break;
                case EventType::arrival :
                    p.calculate_likelihood (e, _arrival_error);
                    break;
                case EventType::departure :
                    p.calculate_likelihood (e, _departure_error);
                    break;
            }

            std::cout << " => l(Y|Xi) = " << exp (p.get_ll ());
        }

        bad_sample = true;
        _Neff = 0;
        double sumlh = std::accumulate (_state.begin (), _state.end (), 0.0,
                                        [](double a, Particle& p) {
                                            return a + p.get_weight () * exp (p.get_ll ());
                                        });

        std::cout << "\n   -> sum(l(y|x)) = " << sumlh;
        // if no likelihoods are that big, give up
        if (sumlh < 1e-10) return;

        // normalize weights
        std::vector<double> wt;
        wt.reserve (_state.size ());
        for (auto& p : _state)
        {
            p.set_weight (p.get_weight () * exp (p.get_ll ()) / sumlh);
            wt.push_back (p.get_weight ());
        }
        double sumwt = std::accumulate (wt.begin (), wt.end (), 0.0);

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
        if (_Neff >= (_N / 4)) return;

        std::cout << " -> resampling";
        select (rng);

    }

    void Vehicle::mutate2 (RNG& rng)
    {
//         if (_newtrip)
//         {
//             initialize (rng);
//             return;
//         }

//         if (_complete || !valid () || _delta == 0) return;

//         // There probably need to be a bunch of checks here ...

//         // std::cout << "\n v" << _vehicle_id 
//         //     << " - M = " << _trip->stops ().size ()
//         //     << "; L = " << _trip->shape ()->segments ().size ();
        
//         // do the transition ("mutation")
//         int ncomplete = 0;
//         // std::cout << " -> travel ... \n";
//         for (auto p = _state.begin (); p != _state.end (); ++p)
//         {
//             if (p->is_complete ())
//             {
//                 ncomplete++;
//             }
//             else
//             {
//                 // std::cout << "\r";
//                 p->travel (_delta, rng);
//             }
//         }
//         // std::cout << "ok";


//         if (ncomplete == _state.size ())
//         {
//             _complete = true;
//             return;
//         }

// // #if WRITE_PARTICLES
// //         std::vector<ShapePt>* path = &(_trip->shape ()->path ());
// //         for (auto p = _state.begin (); p != _state.end (); ++p)
// //         {
// //             p->calculate_likelihood (_position, path, _gpserror);
// //             std::ostringstream fname;
// //             fname << "history/vehicle_" << _vehicle_id << "_proposals.csv";
// //             std::ofstream fout;
// //             fout.open (fname.str ().c_str (), std::ofstream::app);
// //             double d (p->get_distance ());
// //             latlng ppos (_trip->shape ()->coordinates_of (d));
// //             fout << _timestamp << ","
// //                 << _trip->trip_id () << ","
// //                 << p->get_distance () << ","
// //                 << p->get_speed () << ","
// //                 << p->get_acceleration () << ","
// //                 << std::setprecision(15)
// //                 << p->get_ll () << ","
// //                 << ppos.latitude << "," << ppos.longitude << "\n";
// //             fout.close ();
// //         }
// // #endif

// #if SIMULATION
//         // PRIOR model eval stuff
//         double prior_mse = 0.0; // = SUM [W * d(h(X), Y)^2]
//         double dxy;
//         latlng hx;
//         for (auto p = _state.begin (); p != _state.end (); ++p)
//         {
//             dxy = p->get_distance ();
//             hx = _trip->shape ()->coordinates_of (dxy);
//             prior_mse += p->get_weight () * 
//                 pow(distanceEarth (_position, hx), 2);
//         }

//         // prior speed variance
//         double prior_speed = speed ();
//         double prior_speed_var;
//         prior_speed_var = std::accumulate (_state.begin (), _state.end (), 0.0, 
//                                           [&prior_speed](double a, Particle& p) {
//                                             return a + p.get_weight () * pow(p.get_speed () - prior_speed, 2);
//                                           });
// #endif
//         // update
//         // std::cout << " -> select ... ";
//         select (rng);
//         // std::cout << "ok";

//         // (re)initialize if the particle sample is bad
//         bool isbad = bad_sample;
//         if (bad_sample)
//         {
//             initialize (rng);
//         }
//         else
//         {
//             // NETWORK STUFF
//             double dmin = _trip->shape ()->path ().back ().distance;
//             for (auto p = _state.begin (); p != _state.end (); ++p)
//             {
//                 if (p->get_distance () < dmin)
//                 {
//                     dmin = p->get_distance ();
//                 }
//             }
//             std::vector<ShapeSegment>& segs = _trip->shape ()->segments ();
//             int m = find_stop_index (dmin, &(_trip->stops ()));
//             int l = find_segment_index (dmin, &segs);


//             // update segment travel times for intermediate ones ...
//             double tt, ttp, err;
//             int n;
//             while (_current_segment < m)
//             {
//                 // get the average travel time for particles along that segment
//                 tt = 0.0;
//                 n = 0;
//                 for (auto p = _state.begin (); p != _state.end (); ++p)
//                 {
//                     ttp = p->get_travel_time (_current_segment) * p->get_weight ();
//                     if (ttp > 0)
//                     {
//                         tt += ttp;
//                         n++;
//                     }
//                 }
//                 if (n < _N)
//                 {
//                     _current_segment++;
//                     continue;
//                 }

//                 err = std::accumulate (_state.begin (), _state.end (), 0.0,
//                                        [=](double a, Particle& p) {
//                                             return a + p.get_weight () * pow(p.get_travel_time (_current_segment) - tt, 2);
//                                        });

//                 // if the error is effectively 0 ...
//                 if (err < 0.001) err = 10.0;

//                 _segment_travel_times.at (_current_segment) = round (tt);
//                 segs.at (_current_segment).segment->push_data (tt, err, _timestamp);
//                 _current_segment++;
//             }
            
//             // NOTE: need to ignore segment if previous segment travel time is 0
//             // (i.e., can't be sure that the current segment travel time is complete)
            
            
//         }
//         // std::cout << " -> nw";

// #if SIMULATION
//         // POSTERIOR model eval stuff
//         double posterior_mse = 0.0; // = SUM [W * d(h(X), Y)^2]
//         for (auto p = _state.begin (); p != _state.end (); ++p)
//         {
//             dxy = p->get_distance ();
//             hx = _trip->shape ()->coordinates_of (dxy);
//             posterior_mse += p->get_weight () * 
//                 pow(distanceEarth (_position, hx), 2);
//         }

//         // posterior speed variance
//         double post_speed = speed ();
//         double post_speed_var;
//         post_speed_var = std::accumulate (_state.begin (), _state.end (), 0.0, 
//                                           [&post_speed](double a, Particle& p) {
//                                             return a + p.get_weight () * pow(p.get_speed () - post_speed, 2);
//                                           });

//         std::ostringstream mename;
//         mename << "modeleval/vehicle_" << _vehicle_id << ".csv";
//         std::ofstream modeleval;
//         modeleval.open (mename.str ().c_str (), std::ofstream::app);
//         double atd = _trip->shape ()->distance_of (_position);
//         latlng closest_pt = _trip->shape ()->coordinates_of (atd);
//         double ctd = distanceEarth (_position, closest_pt);
//         double sumwt = std::accumulate (_state.begin (), _state.end (), 0.0, [](double a, Particle& p) {
//             return a + p.get_weight ();
//         });
//         double meanwt = sumwt / (double) _N;
//         double varwt = std::accumulate (_state.begin (), _state.end (), 0.0, [&meanwt](double a, Particle& p) {
//             return a + pow (p.get_weight () - meanwt, 2);
//         });
//         varwt /= _N;
//         modeleval << _vehicle_id 
//             << "," << _trip->trip_id ()
//             << "," << _timestamp
//             << "," << prior_mse 
//             << "," << posterior_mse
//             << "," << post_speed
//             << "," << prior_speed_var
//             << "," << post_speed_var
//             << "," << ctd
//             << "," << _Neff
//             << "," << (resample ? 1 : 0)
//             << "," << resample_count
//             << "," << isbad
//             << "\n";
//         modeleval.close ();
// #endif


// // #if WRITE_PARTICLES
// //         for (auto p = _state.begin (); p != _state.end (); ++p)
// //         {
// //             p->calculate_likelihood (_position, path, _gpserror);
// //             std::ostringstream fname;
// //             fname << "history/vehicle_" << _vehicle_id << ".csv";
// //             std::ofstream fout;
// //             fout.open (fname.str ().c_str (), std::ofstream::app);
// //             double d (p->get_distance ());
// //             latlng ppos (_trip->shape ()->coordinates_of (d));
// //             fout << _timestamp << ","
// //                 << _trip->trip_id () << ","
// //                 << std::setprecision(15)
// //                 << _position.latitude << "," << _position.longitude << ","
// //                 << std::setprecision(6)
// //                 << p->get_distance () << ","
// //                 << p->get_speed () << ","
// //                 << p->get_acceleration () << ","
// //                 << std::setprecision(15)
// //                 << p->get_ll () << ","
// //                 << ppos.latitude << "," << ppos.longitude << "\n";
// //             fout.close ();

// //             std::ostringstream fname2;
// //             fname2 << "history/vehicle_" << _vehicle_id << "_particles.csv";
// //             fout.open (fname2.str ().c_str (), std::ofstream::app);
// //             fout << _timestamp;
// //             for (auto ati = p->get_arrival_times ().begin (); ati != p->get_arrival_times ().end (); ++ati)
// //             {
// //                 fout << "," << *ati;
// //             }
// //             fout << "\n";
// //             fout.close ();
// //         }
// // #endif


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



    void Particle::travel (unsigned delta, RNG& rng)
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
        unsigned int m (find_stop_index (distance, stops));
        // std::cout << " [M=" << M << ",m=" << m << "], ";
        if (m == M-1) 
        {
            distance = Dmax;
            complete = true;
            return;
        }

        double next_stop_d = stops->at (m+1).distance;
        
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
        if (distance == stops->at (m).distance)
        {
            if (rng.runif () < 0.05)
            {
                double w = vehicle->gamma () - vehicle->dwell_time () * log (rng.runif ());
                delta = fmax (0, delta - round (w));
                // we don't want this to affect the speed
            }
        }
        // else if (distance == segments->at (l).distance &&
        //          rng.runif () < 0.05)
        // {
        //     double w = - log (rng.runif ()) * delta;
        //     delta = fmax (0, delta - round (w));
        //     if (tt.at (l) >= 0)
        //         tt.at (l) = tt.at (l) + 1;
        // }
        else if (rng.runif () < 0.01)
        {
            // a very small chance for particles to remain stationary
            if (tt.at (l) >= 0)
                tt.at (l) = tt.at (l) + delta;
            delta = 0;

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
        double vmax = rng.runif () < 0.05 ? 30.0 : 15.0;
        while (distance < Dmax && delta > 0.0)
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
            delta--;
            if (tt.at (l) >= 0)
                tt.at (l) = tt.at (l) + 1;
            
            if (l < L-1 && distance >= next_segment_d)
            {
                // reaching intersection ... 
                l++;
                tt.at (l) = 0;
                next_segment_d = (l+1 >= L-1) ? Dmax : segments->at (l+1).distance;
                vmax = rng.runif () < 0.05 ? 30.0 : 15.0;
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
                // about to reach a stop ... slow? stop? just drive past?
                m++; // the stop we are about to reach
                at.at (m) = vehicle->timestamp () - delta;
                if (m >= M-1)
                {
                    distance = next_stop_d;
                    break;
                }
                if (rng.runif () < vehicle->pr_stop ())
                {
                    // stop dwell time ~ Exp(tau = 10)
                    double dwell = vehicle->gamma () - vehicle->dwell_time () * log (rng.runif ());
                    delta = fmax(0, delta - dwell);
                    distance = next_stop_d;
                    // speed = 0.0;
                    // accelerating = 5.0 + rng.runif () * 10.0;
                    // acceleration = 2.0 + rng.rnorm () * vehicle->system_noise ();
                }
                next_stop_d = stops->at (m+1).distance;
                continue;
            }
        }

        // almost at next stop ...
        if (m < M-1 &&
            next_stop_d - distance < 20)
        {
            distance = next_stop_d;
            m++;
            at.at (m) = vehicle->timestamp ();
            if (m < M-1)
            {
                next_stop_d = stops->at (m+1).distance;
            }
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

        std::cout << " => d(h(x), y) = " << exp (ld) << "m";
    }

    void Particle::calculate_likelihood (Event& e, double error)
    {
        uint64_t t = (e.type == EventType::arrival) ? 
            get_arrival_time (e.stop_index) :
            get_departure_time (e.stop_index);

        if (t == 0) 
        {
            std::cout << " => hasn't "
                << (e.type == EventType::arrival ? "arrived" : "departed");
            // log_likelihood = 0.0;
            return;
        }

        std::cout << " => "
            << (e.type == EventType::arrival ? "arrived" : "departed")
            << " at " << t
            << " (diff: " << (e.timestamp - t) << "s)";

        log_likelihood = - 0.5 * log (2 * M_PI) - log (error) - 
            pow (e.timestamp - t, 2) / 2 / pow (error, 2);
        
    }

    void Particle::calculate_arrival_likelihood (int index, uint64_t time, double error)
    {
        // particle's arrival time at stop INDEX
        if (get_arrival_time (index) == 0) log_likelihood = 0;
        log_likelihood = - 0.5 * log (2 * M_PI) - log (error) - pow (get_arrival_time (index) - time, 2) / 2 / pow (error, 2);
    }
    void Particle::calculate_departure_likelihood (int index, uint64_t time, double error)
    {
        log_likelihood = 0;
        // return - 0.5 * log (2 * pi) - log (error) - pow ()
    }

    void Particle::set_weight (double w)
    {
        weight = w;
    }

}; // namespace Gtfs
