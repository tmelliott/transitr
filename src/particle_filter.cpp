/**
 * Define the particle filter functions.
 */
#include "particle_filter.h"
#include "timing.h"

namespace Gtfs {

    void Vehicle::initialize (Event& e, RNG& rng)
    {
#if VERBOSE > 0
        std::cout << "\n    -> initializing";
#endif

        initialize (rng);
        _timestamp = e.timestamp;
        _delta = 0;
        _Neff = _N;

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
            // _current_delay = e.delay;
        }

        // Trip start time:
        Time tstart (_trip->start_time ());
        Time curtime (_timestamp);
        int time_to_start (tstart - curtime);
        // std::cout << "\n -> " << time_to_start << "s until trip is scheduled to begin\n";

        double d, u;
        int si;
        Particle* p;
        auto segs = _trip->shape ()->segments();
        for (int i=0; i<_N; ++i)
        {
            if (e.type == EventType::gps)
            {
                if (time_to_start > 60)
                {
                    d = 0.;
                }
                else
                {
                    // initialize each particle within 100m of obs
                    u = rng.runif ();
                    d = (dist < 0 ? u * dmax : fmax(0, fmin(dmax, u * 200 - 100 + dist)));
                }

                si = find_segment_index (d, &segs);
                _state.emplace_back (
                    d,
                    segs.at (si).segment->sample_speed (rng),
                    rng.rnorm () * _systemnoise, // acceleration, not used (currently)
                    this
                );
            }
            else
            {
                // stick the bus AT the stop
                si = find_segment_index (dist, &segs);
                _state.emplace_back (
                    dist,
                    segs.at (si).segment->sample_speed (rng), // speed
                    rng.rnorm () * _systemnoise, // acceleration, not used (currently)
                    this
                );

                // point to the particle
                p = &(_state.back ());
                p->bus_stop (e.timestamp, rng);
                if (e.type == EventType::departure)
                {
                    // shift at <- dt
                    int dwell (p->get_departure_time (si) - p->get_arrival_time (si));
                    p->set_arrival_time (si, p->get_arrival_time (si) - dwell);
                    p->set_departure_time (si, p->get_departure_time (si) - dwell);
                }

                // and set the travel time for the segment to zero
                p->init_travel_time (si);
            }
        }

        double initwt = 1.0 / (double) _N;
        for (auto p = _state.begin (); p != _state.end (); ++p) p->set_weight (initwt);

#if SIMULATION
        this->store_state ("initialize");
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
        _segment_speed_avgs.clear ();
        _stop_arrival_times.clear ();
        _stop_departure_times.clear ();
        // _tt_state.clear ();
        // _tt_cov.clear ();

        _segment_speed_avgs.resize (_trip->shape ()->segments ().size (), 0);
        _stop_arrival_times.resize (_trip->stops ().size (), 0);
        _stop_departure_times.resize (_trip->stops ().size (), 0);

        // alright: initialize these using the schedule
        // _tt_state.resize (_trip->stops ().size (), 0.0);
        // _tt_cov.resize (_trip->stops ().size (),
        //                 std::vector<double> (_trip->stops ().size (), 0.0));
        // for (int i=0; i<_trip->stops ().size (); i++)
        // {
        //     _tt_cov.at (i).at (i) = (i+1) * 30; // 5 min error + 30 seconds per stop
        // }
    }

    void Vehicle::mutate (RNG& rng, Gtfs* gtfs)
    {
        int nnew = time_events.size () - current_event_index;
#if VERBOSE > 0
        std::cout << "\n\n+ vehicle " << _vehicle_id << ": "
            << nnew << " new events";
#endif
        if (nnew == 0) return;

#if VERBOSE > 0
        if (current_event_index == 0 || _timestamp == 0)
        {
            std::cout << "\n    [uninitialized]";
        }
        else
        {
            if (current_event_index > 1)
            {
                for (int i=0; i<(current_event_index-1); ++i)
                {
                    std::cout << "\n    [" << Time (time_events.at (i).timestamp) << "] ";
                    time_events.at (i).print ();
                }
            }
            std::cout << "\n    [" << Time (_timestamp) << "] ";
            if (_trip == nullptr)
            {
                std::cout << "the trip is null ... ";
            }
            else
            {
                std::cout
                    << _trip->route ()->route_short_name ()
                    << " (" << _trip->stops ().at (0).departure_time << "): ";
            }
            time_events.at (current_event_index - 1).print ();
            std::cout << "\n ------------------";
        }
#endif


        // repeat until there are no more events
        while (current_event_index < time_events.size ())
        {
            // std::cout << std::endl << " ++++ there are " << time_events.size ()
            //     << " events; requesting index " << current_event_index << std::endl;
            auto e = time_events.at (current_event_index);
            _latest_event = &(time_events.at (current_event_index));

            if (_trip == nullptr || _trip->trip_id () != e.trip_id)
            {
                // unload old trip
                if (_trip != nullptr) _trip->unload (true);

                // assign trip <--> vehicle
                try
                {
                    set_trip (gtfs->find_trip (e.trip_id));
                }
                catch (int e)
                {
                    std::cout << "\n - Another vehicle is already assigned to that trip! Skipping event ...";
                    current_event_index++;
                    _trip = nullptr;
                    continue;
                }
                _newtrip = _trip != nullptr;

                _previous_state.clear ();
                _previous_ts = 0;
                _timestamp = 0;
                _state.clear ();
                estimated_dist = 0.0;
                _current_delay = 0;

                // reset events!!
                new_events.clear ();
                // move future events (including current) into new_events
                for (int i=current_event_index; i<time_events.size (); i++) new_events.push_back (std::move (time_events.at (i)));
                time_events.clear ();
                // move back into events
                for (int i=0; i<new_events.size (); i++) time_events.push_back (std::move (new_events.at (i)));
                new_events.clear ();
                current_event_index = 0;
            }

            if (_trip == nullptr)
            {
                throw std::runtime_error ("Trip not found");
            }

#if VERBOSE > 0
            std::cout << "\n    [" << e.timestamp << "] "
                << _trip->route ()->route_short_name ();
            if (_trip->stops ().size () == 0)
            {
                std::cout << " -> has no stops ... :/";
            }
            else
            {
                std::cout
                    << " (" << _trip->stops ().at (0).departure_time
                    << ", " << _trip->stops ().size () << " stops"
                    << "): ";
            }
            e.print ();
#endif

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
#if VERBOSE > 0
                        std::cout << " -> distance to route = " << dist_to_route << "m";
#endif
                        if (dist_to_route > 50)
                        {
                            current_event_index++;
                            return;
                        }

                        // if the estimated distance is less than the previous
                        // estimated distance, do something about that
                        auto pos = latlng (e.position.latitude, e.position.longitude);
                        auto est_dist = _trip->shape ()->distance_of (pos);
                        // std::cout << "\n-> from " << estimated_dist << " to " << est_dist;
                        if (est_dist - estimated_dist < -5.0)
                        {
                            // looks like we've gone backwards ...
                            // ... uh oh
                            // 1. skip this observation
                            // 2. undo the entire state, and repredict state using newest obs
                            // 3. revert particle weights and recalculate based on this new observation
                            //    (this is probably the easiest, tbh)
                            this->revert_state ();
                        }
                        break;
                    }
                case EventType::arrival :
                    {
                        // this is tricky ...
                        _trip->set_arrival_time (e.stop_index, e.timestamp);
                        // if (e.stop_index == _trip->stops ().size () - 1)
                        // {
                        //     _complete = true;
                        // }
                        break;
                    }
                case EventType::departure :
                    {
                        // no checks
                        _trip->set_departure_time (e.stop_index, e.timestamp);
                        // if (e.stop_index == _trip->stops ().size () - 1)
                        // {
                        //     _complete = true;
                        // }
                        break;
                    }
            }


            // pass the event data
            if (e.type == EventType::gps)
            {
                _position = latlng (e.position.latitude, e.position.longitude);
                estimated_dist = _trip->shape ()->distance_of (_position);
            }
            else
            {
                _stop_index = e.stop_index;
                _current_delay = e.delay;
                // ----------------------------------- CHECK STOP INDEX
#if VERBOSE > 0
                std::cout << " - there are " << _trip->stops ().size ()
                    << " stops " << " (d = "
                    << _trip->stops ().at (_stop_index).distance << "m)";
#endif
            }
            // std::cout << "\n - the vehicle's current delay is " << _current_delay << "s";

            if (!_skip_observation)
            {
                _previous_state = _state;
                mutate_to (e, rng);
            }
            if (_skip_observation)
            {
                std::cout << "\n - skipping observation ...";
            }
            _skip_observation = false;

            // if the current iteration fails, start again from here
            if (bad_sample) initialize (e, rng);

            {
#if VERBOSE > 0
                std::cout << "\n    ** estimating speeds";
#endif
                // NETWORK STUFF
                std::vector<ShapeSegment>& segs = _trip->shape ()->segments ();

                // update segment travel times for intermediate ones ...
                double sp, ttp, err;
                int n;
                int M = segs.size ();
                double segmin, segmax;

#if SIMULATION
                // std::ofstream fout;
                // fout.open ("particle_travel_times.csv", std::ofstream::app);
#endif

                double L;
                for (_current_segment=0; _current_segment<M; _current_segment++)
                {
                    if (_segment_speed_avgs.at (_current_segment) > 0) continue;
                    L = segs.at (_current_segment).segment->length ();
                    segmin = L / segs.at (_current_segment).segment->max_speed ();
                    segmax = L * 5.;

                    // get the average speed for particles along that segment
                    sp = 0.0;
                    n = 0;
                    for (auto p = _state.begin (); p != _state.end (); ++p)
                    {
                        if (p->get_travel_time (_current_segment) < segmin) continue;
                        if (p->get_travel_time (_current_segment) > segmax) continue;
                        if (p->get_segment_index () == _current_segment) continue;
                        sp += p->get_weight () *
                            L / p->get_travel_time (_current_segment);
                        n++;
                    }
                    if (n < _N ||
                        sp < 0.5 ||
                        sp > segs.at (_current_segment).segment->max_speed ())
                    {
                        continue;
                    }

#if VERBOSE > 0
                    std::cout << "\n  - segment "
                        << (_current_segment + 1)
                        << " of " << _segment_speed_avgs.size ();
#endif

                    // this is the variance
                    err = std::accumulate (
                        _state.begin (),
                        _state.end (),
                        0.0,
                        [=](double a, Particle& p) {
                            return a + p.get_weight () *
                                pow (L / p.get_travel_time (_current_segment) - sp, 2);
                        }
                    );

                    _segment_speed_avgs.at (_current_segment) = (sp);

                    // if (_stop_arrival_times.at (_current_segment+1) > 0 &&
                    //     _stop_departure_times.at (_current_segment) > 0)
                    // {
                    //     sp = _stop_arrival_times.at (_current_segment+1) -
                    //         _stop_departure_times.at (_current_segment);
                    //     err = 3.;
                    // }

                    segs.at (_current_segment).segment->push_data (sp, err, _timestamp);
#if VERBOSE > 0
                    std::cout << ": " <<  (sp) << " (" << err << ")";
#endif
#if SIMULATION
                    // for (auto p = _state.begin (); p != _state.end (); ++p)
                    // {
                    //     fout << _timestamp
                    //         << "," << _vehicle_id
                    //         << "," << segs.at (_current_segment).segment->segment_id ()
                    //         << "," << p->get_travel_time (_current_segment)
                    //         << "," << p->get_weight ()
                    //         << "\n";
                    // }
#endif
                }

                // NOTE: need to ignore segment if previous segment travel time is 0
                // (i.e., can't be sure that the current segment travel time is complete)


#if SIMULATION
                // fout.close ();
#endif
            }

            current_event_index++;
        }

        // is the vehicle finished?
        if (_complete)
        {
            // delete its particles, because we don't need them
            _state.clear ();
        }



        return;
    }

    void Vehicle::mutate_to (Event& e, RNG& rng)
    {
        bool log = false;
#if VERBOSE > 0
        log = true;
#endif
        if (log) std::cout << "\n -> mutating to event ...";

        if (_newtrip || bad_sample)
        {
            if (log) std::cout << " initializing";
            initialize (e, rng);
            return;
        }
        if (log) std::cout << "here we go!";
        bad_sample = true;

        _delta = e.timestamp - _timestamp;
        _timestamp = e.timestamp;

        if (e.type == EventType::arrival)
        {
            _stop_arrival_times.at (e.stop_index) = e.timestamp;
        }
        if (e.type == EventType::departure)
        {
            _stop_departure_times.at (e.stop_index) = e.timestamp;
        }

        // move the particles
        if (_complete || !valid () || _delta == 0)
        {
            if (log)
            {
                std::cout << "\n ... uh oh ... ";
                if (_complete) std::cout << "complete!";
                if (!valid ()) std::cout << "invalid!";
                if (_delta == 0) std::cout << "delta = 0!";
            }
            return;
        }
#if VERBOSE > 0
        std::cout << "\n     + " << _delta << " seconds";
#endif

        bool all_complete = true;
        double dbar = 0.0, vbar = 0.0,
            dbar2 = 0.0, vbar2 = 0.0, ddbar = 0.0,
            dtbar = 0;
        int firstN = 20;
        int pi = 0;
        for (auto& p : _state)
        {
            pi++;
#if VERBOSE > 0
            if (firstN > 0)
                std::cout << "\n      ["
                    << std::setw (5) << round (p.get_distance ())
                    << ", " << std::setw (4) << (round (p.get_speed () * 10) / 10)
                    << ", " << std::setw (2) << (p.get_stop_index () + 1)
                    << ", " << std::setw (2) << (p.get_segment_index () + 1)
                    << "]";
            // else
            //     std::cout << "\r                                                                               "
            //         << "\r transitioning particle ... " << pi;
#endif
            dbar += p.get_distance () * p.get_weight ();
            vbar += p.get_speed () * p.get_weight ();

            p.travel (_delta, e, rng);
#if VERBOSE > 0
            // if (firstN == 0)
            //     std::cout << " complete ...";
#endif


            // if any aren't complete, prevent vehicle from finishing trip
            all_complete = false;
#if VERBOSE > 0
            if (firstN > 0) {
                std::cout << " -> ["
                    << std::setw (5) << round (p.get_distance ())
                    << ", " << std::setw (4) << (round (p.get_speed () * 10) / 10)
                    << ", " << std::setw (2) << (p.get_stop_index () + 1)
                    << ", " << std::setw (2) << (p.get_segment_index () + 1)
                    << "]";
            }
#endif

            dbar2 += p.get_distance () * p.get_weight ();
            vbar2 += p.get_speed () * p.get_weight ();

            if (p.is_complete ())
            {
#if VERBOSE > 0
                if (firstN > 0) {
                    // firstN--;
                    std::cout << " -> complete";
                }
#endif
                // continue;
            }

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
#if VERBOSE > 0
                        if (firstN > 0) {
                            std::cout << " -> d(Y,Xi) = "
                                << std::setw (5) << round (distanceEarth (e.position, pos));
                        }
#endif
                        break;
                    }
                case EventType::arrival :
                    p.calculate_likelihood (e, _arrival_error);
                    dtbar += (int)(p.get_arrival_time (e.stop_index) - e.timestamp) * p.get_weight ();
#if VERBOSE > 0
                        if (firstN > 0) {
                            int art = p.get_arrival_time (e.stop_index);
                            std::cout << ", arr = "
                                << Time (art)
                                << " -> d(A,Ai) = "
                                << std::setw (3) << (int)(p.get_arrival_time (e.stop_index) - e.timestamp);
                        }
#endif
                    break;
                case EventType::departure :
                    p.calculate_likelihood (e, _departure_error);
                    dtbar += (int)(p.get_departure_time (e.stop_index) - e.timestamp) * p.get_weight ();
#if VERBOSE > 0
                        if (firstN > 0) {
                            int art = p.get_arrival_time (e.stop_index);
                            std::cout << ", dep = "
                                << Time (art)
                                << " -> d(D,Di) = "
                                << std::setw (3) << (int)(p.get_departure_time (e.stop_index) - e.timestamp);
                        }
#endif
                    break;
            }

#if VERBOSE > 0
            if (firstN > 0) {
                firstN--;
                std::cout << " => l(Y|Xi) = " << exp (p.get_ll ());
            } else {
                // std::cout << "likelihood calculated ... done!";
            }
#endif
        }

        latlng px = latlng ();
#if VERBOSE > 0
        std::cout <<
            "\n    ========================================================================="
            << std::endl
            << "      [" << dbar << ", " << vbar << "] -> "
            << "[" << dbar2 << ", " << vbar2<< "] => ";

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
#endif

        if (all_complete) _complete = true;

        // _Neff = 0;
        double sumlh = std::accumulate (_state.begin (), _state.end (), 0.0,
                                        [](double a, Particle& p) {
                                            return a + p.get_weight () * exp (p.get_ll ());
                                        });

#if SIMULATION
        this->store_state ("mutate");
#endif

#if VERBOSE > 0
        std::cout << "\n   -> sum(l(y|x)) = " << sumlh;
#endif
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

#if VERBOSE > 0
        std::cout << "\n         => Posterior: ["
            << dbar << ", " << vbar << "] => ";
        switch (e.type)
        {
            case EventType::gps :
                {
                    std::cout << "(getpx_len="
                        << _trip->shape ()->path ().size () << ")";
                    px = _trip->shape ()->coordinates_of (dbar);
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
#endif

        if (fabs (1 - sumwt) > 1e-4) return;
        bad_sample = false;

        // calcualte Neff.
        double sumwt2 = std::accumulate(wt.begin (), wt.end (), 0.0,
                                        [](double a, double b) {
                                            return a + pow(b, 2);
                                        });

        _Neff = pow (sumwt2, -1);
#if VERBOSE > 0
        std::cout << "\n   -> Neff = " << _Neff;
#endif

        if (_Neff < (_N / 4)) action = "resample";
#if SIMULATION
        this->store_state ("update");
#endif


        if (_Neff >= (_N / 4) || _complete) return;

#if VERBOSE > 0
        std::cout << " -> resampling";
#endif
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
#if VERBOSE > 3
        std::cout << "\n ===> Transition particle\n"
            << " * Initial state: [d=" << distance
            << ", s=" << speed
            << ", j=" << stop_index
            << ",l=" << segment_index
            << "]";

        std::cout << "\n * Delta = " << delta;

        std::cout << "\n * Event = ";
        if (e.type == EventType::gps)
        {
            std::cout << "GPS observation";
        }
        else
        {
            std::cout << (e.type == EventType::arrival ? "arrival" : "departure")
                << " at stop " << (e.stop_index + 1);
        }
#endif
        // Shape is:
        Shape* shape = vehicle->trip ()->shape ();
        std::vector<ShapeNode>& nodes = shape->nodes ();
        std::vector<StopTime>& stops = vehicle->trip ()->stops ();
        std::vector<ShapeSegment>& segments = shape->segments ();

        int M, L;
        M = stops.size ();
        L = nodes.size ();
#if VERBOSE > 3
        std::cout << "\n * Shape: " << shape->path ().size ()
            << " with " << L << " nodes and " << M << " stops";
#endif

        if (stop_index == M - 1 || segment_index == L - 1)
        {
            complete = true;
        }

        if (delta == 0 || complete || !vehicle || !vehicle->trip () ||
            !vehicle->trip ()->shape ())
        {
#if VERBOSE > 3
            std::cout << "\n ** nothing to mutate, skipping particle\n";
#endif
            return;
        }

        // Necessary to clear any extrapolated stops
        at.at (stop_index + 1) = 0;
        dt.at (stop_index + 1) = 0;

        double Dmax (shape->path ().back ().distance);
#if VERBOSE > 3
        std::cout << "\n * Route distance = " << Dmax << "m";
#endif
        if (distance >= Dmax)
        {
            distance = Dmax;
            complete = true;
#if VERBOSE > 3
            std::cout << " -> particle has finished\n\n";
#endif
            return;
        }

        ShapeNode* next_node;
        next_node = &(nodes.at (segment_index + 1));
#if VERBOSE > 3
        std::cout << "\n\n * Next node: "
            << (segment_index + 1)
            << " (" << next_node->distance << "m)";
#endif

        StopTime* next_stop;
        next_stop = &(stops.at (stop_index + 1));
#if VERBOSE > 3
        std::cout << "\n * Next stop: "
            << (stop_index + 1)
            << " (" << next_stop->distance << "m)";
#endif

        if (next_node->node->node_type () == 0)
        {
            if (next_node->node->node_id () == next_stop->stop->node ()->node_id ())
            {
#if VERBOSE > 3
                std::cout << "\n   -> next node is the next stop!";
#endif
            }
            else
            {
                // something odd happening and stop/node aren't matching up :(
#if VERBOSE > 3
                std::cout << "\n   -> next node is a stop, but not the right one!!!\n\n";
#endif
                return;
            }
        }

        // if at a stop, we must wait there while people get on/off
        if (distance == nodes.at (segment_index).distance)
        {
#if VERBOSE > 3
            std::cout << "\n   -> currently still at node " << segment_index;
#endif
            double wait (0);
            if (nodes.at (segment_index).node->node_type () == 0)
            {
#if VERBOSE > 3
                std::cout << " which is a bus stop";
#endif
                if (segment_index == 0)
                {
                    Time& sched = stops.at (0).departure_time;
                    double avgdelay = stops.at (0).average_delay;
                    double sddelay = stops.at (0).sd_delay;

                    int depi;
                    if (sddelay == 0)
                    {
                        depi = rng.runif () * 90 - 30; // [-30, 60] second delay
                    }
                    else
                    {
                        depi = round (rng.rnorm () * sddelay + avgdelay);
                    }

                    if (Time (e.timestamp) > sched + depi)
                    {
                        wait = 0;
                    }
                    else
                    {
                        wait = (sched + depi) - Time (e.timestamp);
                    }
                }
                else if (stops.at (stop_index).departure_time > stops.at (stop_index).arrival_time)
                {
                    // iteration start time = e.timestamp - delta
                    // scheduled departure = stops.at.dep time
                    wait = fmax (
                        0.0,
                        stops.at (stop_index).departure_time -
                            (Time (e.timestamp - delta)) +
                            rng.rnorm () * 5.0
                    );
                }
                else if (rng.runif () < vehicle->params ()->pr_stop)
                {
                    wait = vehicle->gamma () +
                        vehicle->dwell_time () +
                        vehicle->dwell_time_var () * rng.rnorm ();
                }
                else
                {
                    wait = 0.0;
                }
            }
            else
            {
#if VERBOSE > 3
                std::cout << " which is a generic intersection";
#endif
                if (rng.runif () < vehicle->params ()->pr_stop)
                {
                    // wait = - vehicle->dwell_time () * log (rng.runif ());
                    // ignore segment intersection times for now ...
                    wait = 0.0;
                }
                else
                {
                    wait = 0.0;
                }
            }
#if VERBOSE > 3
            std::cout << " -> waiting " << round (wait) << "s";
#endif
            delta -= round (wait);
        }

        // if it's the last stop, gotta go to the end
        bool complete_route
            (e.type != EventType::gps && e.stop_index == stops.size () - 1);

        if (delta < 0 && !complete_route) return;

        if (vehicle->params ()->noise_model == 0)
        {
#if VERBOSE > 3
            std::cout << "\n\n * Using noise model 0, adjusting speed:"
                << "\n   " << speed << " -> ";
#endif
            // add noise once per iteration,
            // or when passing segment/stop
            double vel = speed + rng.rnorm () * vehicle->system_noise ();
            int nn = 100;
            double vmax = segments.at (segment_index).segment->max_speed ();
            while (vel <= 0 || vel > vmax)
            {
                vel = speed + rng.rnorm () * vehicle->system_noise ();
                if (nn-- == 0) vel = rng.runif () * vmax;
            }
            speed = vel;
#if VERBOSE > 3
            std::cout << speed;
#endif
        }

        // now begin travelling
        double node_dist, node_eta;
        while ((complete_route && !complete) || delta > 0)
        {
            // noise model 0
            // jump straight to node/
            if (vehicle->params ()->noise_model == 0)
            {
                node_dist = next_node->distance - distance;
                node_eta = node_dist / speed;
#if VERBOSE > 3
                std::cout << "\n\n * Next node is "
                    << "(" << next_node->distance << " - " << distance << ") = "
                    << node_dist << "m away"
                    << "\n   -> ETA of " << node_eta << "s, delta = "
                    << delta << ", j = " << stop_index << ", l = " << segment_index;
#endif
                // if trying to complete route, don't use this method
                if (node_eta > delta && !complete_route)
                {
                    distance += delta * speed;
                    tt.at (segment_index) += delta;
                    delta = 0;
#if VERBOSE > 3
                    std::cout << "\n   -> stopped at " << distance << "m";
#endif
                }
                else if (segment_index == L - 2)
                {
                    delta -= node_eta;
#if VERBOSE > 3
                    std::cout << " -> arrived at last node";
#endif
                    tt.at (segment_index) += node_eta;
                    distance = Dmax;
                    stop_index++;
                    at.at (stop_index) = vehicle->timestamp () - delta;
#if VERBOSE > 3
                    std::cout << " at " << at.at (stop_index);
                    std::cout << "\n   -> last segment travel time is "
                        << tt.back ()
                        << " (segment " << segment_index << ")";
#endif
                    complete = true;
                    segment_index++;
                    delta = 0;
                    break;
                }
                else
                {
                    delta -= node_eta;
                    tt.at (segment_index) += node_eta;
                    distance = next_node->distance;
                    segment_index++;
                    next_node = &(nodes.at (segment_index + 1));


                    bool is_stop;
                    is_stop = next_node->node->node_type () == 0;
                    if (is_stop)
                    {
                        stop_index++;
                        next_stop = &(stops.at (stop_index + 1));
                    }

                    // adjust speed
                    // if (segments.at (segment_index).segment->travel_time () > 0 &
                    //     segments.at (segment_index).segment->uncertainty () > 0 &
                    //     segments.at (segment_index).segment->uncertainty () < segments.at (segment_index).segment->length () * 2)
                    // {
                    //     speed = segments.at (segment_index).segment->sample_speed (rng);
                    // }
                    // else
                    // {
                    //     // set speed to scheduled speed +- some noise
                    //     if (next_node->node->node_type () == 0 &&
                    //         nodes.at (segment_index).node->node_type () == 0)
                    //     {
                    //         int sched_tt = stops.at (stop_index + 1).arrival_time -
                    //             stops.at (stop_index).departure_time;
                    //         // std::cout << "\n - inter-stop time = " << sched_tt << "s";
                    //         // std::cout << "\n - length = " << segments.at (segment_index).segment->length () << "m";
                    //         speed = segments.at (segment_index).segment->length () / sched_tt;
                    //         // std::cout << "\n - speed = " << speed;
                    //     }
                    // }
#if VERBOSE > 3
                    std::cout << "\n -> setting vehicle speed to " <<
                        speed << " based on segment state";
#endif

#if VERBOSE > 3
                    std::cout << "\n   -> arrival at node " << (segment_index);
#endif
                    tt.at (segment_index) = 0;
                    at.at (stop_index) = vehicle->timestamp () - delta;
#if VERBOSE > 3
                    std::cout << " at " << at.at (stop_index);
#endif

                    // now handle stopping behaviour
                    double dwell = 0.0;
                    if (is_stop)
                    {
                        // if (rng.runif() < vehicle->pr_stop ())

                        // Pr(stop) inverse proportional to speed
                        // [0m/s -> 1, 35m/s -> 0]
                        //if (rng.runif () > vehicle->params ()->pr_stop)
                        if (rng.runif () > 1. - speed / 35.)
                        {
#if VERBOSE > 3
                            std::cout << "\n   -> dwell time at stop: ";
#endif
                            if (stops.at (stop_index).sd_delay > 0)
                            {
                                dwell = 0.0;
                                while (dwell <= 0.0)
                                {
                                    dwell = rng.rnorm () * stops.at (stop_index).sd_delay +
                                        stops.at (stop_index).average_delay;
                                }
                                dwell = fmax (vehicle->gamma (), dwell);
                            }
                            else
                            {
                                dwell = vehicle->gamma () +
                                    vehicle->dwell_time () +
                                    vehicle->dwell_time_var () * rng.rnorm ();
                            }
#if VERBOSE > 3
                            std::cout << round (dwell);
#endif
                        }
                        else
                        {
#if VERBOSE > 3
                            std::cout << "\n   -> skipping stop";
#endif
                        }
                        dt.at (stop_index) = at.at (stop_index) + round (dwell);
#if VERBOSE > 3
                        std::cout << " -> departure time = "
                            << dt.at (stop_index);
#endif
                    }
                    else
                    {
                        if (rng.runif () < vehicle->pr_stop ())
                        {
#if VERBOSE > 3
                            std::cout << "\n   -> queue time at node: ";
#endif
                            dwell = vehicle->dwell_time () +
                                vehicle->dwell_time_var () * rng.rnorm ();
#if VERBOSE > 3
                            std::cout << round (dwell);
#endif
                        }
                        else
                        {
#if VERBOSE > 3
                            std::cout << "\n   -> not stopping";
#endif
                        }
                    }
                    delta -= fmin (delta, round(dwell));
                }

            }
        }

#if VERBOSE > 3
        std::cout << "\n\n * Final state: [d=" << distance
            << ", s=" << speed
            << ", j=" << stop_index
            << ", l=" << segment_index
            << "]";
#endif

        if (e.type != EventType::gps && behind_event (e, delta))
        {
#if VERBOSE > 3
            std::cout << "\n * Event is a stop arrival/departure ... "
                << "extrapolating particle's arrival/departure time";
#endif
            if (at.at (stop_index + 1) == 0)
            {
                double stop_dist, stop_eta;
                stop_dist = fmax (0, next_stop->distance - distance);
                stop_eta = stop_dist / speed;
                at.at (stop_index + 1) = vehicle->timestamp () + stop_eta;
#if VERBOSE > 3
                std::cout << "\n   -> arrival: "
                    << at.at (stop_index + 1);
#endif
            }

            if (e.type == EventType::departure)
            {
                double dwell = 0;
                if (rng.runif () < vehicle->pr_stop ())
                {
                    dwell += vehicle->gamma () +
                        vehicle->dwell_time () +
                        vehicle->dwell_time_var () * rng.rnorm ();
                }
                dt.at (stop_index + 1) = at.at (stop_index + 1) + dwell;
#if VERBOSE > 3
                std::cout << "\n   -> departure: "
                    << dt.at (stop_index + 1);
#endif
            }
        }

#if VERBOSE > 3
        std::cout << "\n\n ==> transition complete.\n"
            << "--------------------------\n\n";
#endif
    }

    // void Particle::old_travel (int delta, Event& e, RNG& rng)
    // {
            // if (complete || !vehicle || !vehicle->trip () ||
            // !vehicle->trip ()->shape ()) return;

    //     // do the particle physics
    //     double Dmax = vehicle->trip ()->shape ()->path ().back ().distance;
    //     if (distance >= Dmax)
    //     {
    //         distance = Dmax;
    //         complete = true;
    //         return;
    //     }

    //     // get STOPS
    //     std::vector<StopTime>* stops;
    //     stops = &(vehicle->trip ()->stops ());
    //     int M (stops->size ());
    //     stop_index = find_stop_index (distance, stops);
    //     // std::cout << " [M=" << M << ",m=" << stop_index << "], ";
    //     if (stop_index == M-1)
    //     {
    //         distance = Dmax;
    //         complete = true;
    //         return;
    //     }

    //     double next_stop_d = stops->at (stop_index + 1).distance;

    //     // get SEGMENTS
    //     std::vector<ShapeSegment>* segments;
    //     segments = &(vehicle->trip ()->shape ()->segments ());
    //     int L (segments->size ());
    //     // std::cout << " [L=" << L;
    //     unsigned int l (segment_index);
    //     // std::cout << ",l=" << l << " --> " << segment_index << "], ";
    //     double next_segment_d;
    //     next_segment_d = (l+1 >= L-1) ? Dmax : segments->at (l+1).distance;

    //     // allow vehicle to remain stationary if at a stop:
    //     // if (distance == stops->at (stop_index).distance)
    //     // {
    //     //     if (rng.runif () < 0.05)
    //     //     {
    //     //         double w = vehicle->gamma () - vehicle->dwell_time () * log (rng.runif ());
    //     //         delta = fmax (0, delta - round (w));
    //     //         // we don't want this to affect the speed
    //     //     }
    //     // }
    //     // else if (distance == segments->at (l).distance &&
    //     //          rng.runif () < 0.05)
    //     // {
    //     //     double w = - log (rng.runif ()) * delta;
    //     //     delta = fmax (0, delta - round (w));
    //     //     if (tt.at (l) >= 0)
    //     //         tt.at (l) = tt.at (l) + 1;
    //     // }


    //     // adjust "delta" if arrival/departure already set
    //     uint64_t tx = vehicle->timestamp () - delta;
    //     if (dt.at (stop_index) > tx)
    //     {
    //         // particle is AT the stop, not yet departed (it's just about to...)
    //         if (e.type == EventType::departure && stop_index == e.stop_index)
    //         {
    //             uint64_t newtime = rng.rnorm () * vehicle->departure_error () + e.timestamp + 0.5;
    //             // departure time has to be >= arrival time
    //             dt.at (stop_index) = fmax (newtime, at.at (stop_index));
    //             delta = (int)(vehicle->timestamp () - dt.at (stop_index));
    //         }
    //         // can't fudge arrival time (in case its based on speed!!!)
    //     }
    //     else if (at.at (stop_index) > tx)
    //     {
    //         // particle has not yet ARRIVED at the stop
    //         if (e.type == EventType::departure && stop_index == e.stop_index)
    //         {
    //             // let's fudge
    //             uint64_t newtime = rng.rnorm () * vehicle->departure_error () + e.timestamp + 0.5;
    //             dt.at (stop_index) = fmax (newtime, at.at (stop_index));
    //             delta = (int)(vehicle->timestamp () - dt.at (stop_index));
    //         }
    //         else if (e.type == EventType::arrival && stop_index == e.stop_index)
    //         {
    //             // resample dwell time (to allow variation)
    //             bus_stop (at.at (stop_index), rng);
    //             delta = (int)(vehicle->timestamp () - dt.at (stop_index));
    //         }
    //         else
    //         {
    //             delta = (int)(vehicle->timestamp () - at.at (stop_index));
    //         }
    //     }

    //     if (distance > stops->at (stop_index).distance && rng.runif () < 0.01)
    //     {
    //         // a very small chance for particles to remain stationary
    //         int tr = rng.runif () * delta + 0.5;
    //         if (tt.at (l) >= 0)
    //             tt.at (l) = tt.at (l) + tr;
    //         delta -= tr;

    //         // speed = 0.0;
    //         // when /not/ at a bus stop, set speed to 0 and wait
    //         // double w = - log (rng.runif ()) * delta;
    //         // delta = fmax (0.0, delta - round (w));
    //         // then the bus needs to accelerate back up to speed ... for how many seconds?
    //         // accelerating = 5.0 + rng.runif () * 10.0;
    //         // acceleration = 2.0 + rng.rnorm () * vehicle->system_noise ();
    //     }


    //     double speed_mean = 10.0;
    //     double speed_sd = 100.0;
    //     // if (segments->at (l).segment->travel_time () > 0 &&
    //     //     segments->at (l).segment->uncertainty () > 0)
    //     // {
    //     //     speed_mean = segments->at (l).segment->get_speed ();
    //     //     speed_sd = - segments->at (l).segment->length () / pow (speed_mean, 2) *
    //     //         segments->at (l).segment->uncertainty ();
    //     //     speed_sd = pow(speed_sd, 0.5);
    //     // }
    //     double vmax = 30; //rng.runif () < 0.5 ? 30.0 : 15.0;

    //     if (vehicle->params ()->noise_model == 0)
    //     {
    //         // add noise once per iteration,
    //         // or when passing segment/stop
    //         double vel = speed + rng.rnorm () * vehicle->system_noise ();
    //         while (vel <= 0 || vel > 30) vel = speed + rng.rnorm () * vehicle->system_noise ();
    //         speed = vel;
    //     }

    //     // while (distance < Dmax && delta > 0.0)
    //     while (behind_event (e, delta))
    //     {
    //         // add system noise to acceleration to ensure speed remains in [0, vmax]
    //         // double accel_prop (-100.0);
    //         // double n = 0;
    //         // if (speed > vmax)
    //         // {
    //         //     speed = rng.runif () * vmax;
    //         // }
    //         // while (speed + accel_prop < 0.0 || speed + accel_prop > vmax && n < 1000)
    //         // {
    //         //     accel_prop = rng.rnorm () * vehicle->system_noise () *
    //         //         (1.0 + (double)n / 100.0);
    //         //     n++;
    //         //     // if (accelerating > 0.0)
    //         //     // {
    //         //     //     accel_prop += acceleration;
    //         //     //     accelerating--;
    //         //     // }
    //         // }

    //         // double v = fmax (0, fmin (30, speed + acceleration));
    //         // double v = speed;
    //         // double vstar = speed + accel_prop;
    //         // double alpha = (pow (v - speed_mean, 2) - pow(vstar - speed_mean, 2)) / (2 * pow (speed_sd, 2));
    //         // alpha = fmin (0, alpha);
    //         // if (rng.runif () < exp (alpha))
    //         // {
    //         //     acceleration = accel_prop;
    //         //     speed = vstar;
    //         // }
    //         // else
    //         // {
    //         //     speed = v;
    //         // }

    //         if (vehicle->params ()->noise_model == 1)
    //         {
    //             double vel = speed + rng.rnorm () * vehicle->system_noise ();
    //             while (vel <= 0 || vel > 30) vel = speed + rng.rnorm () * vehicle->system_noise ();
    //             speed = vel;
    //         }
    //         distance += speed;
    //         if (distance >= Dmax) complete = true;
    //         delta--;

    //         if (tt.at (l) >= 0)
    //             tt.at (l) = tt.at (l) + 1;

    //         if (l < L-1 && distance >= next_segment_d)
    //         {
    //             // reaching intersection ...
    //             l++;
    //             segment_index++;
    //             tt.at (l) = 0;
    //             next_segment_d = (l+1 >= L-1) ? Dmax : segments->at (l+1).distance;
    //             // vmax = rng.runif () < 0.5 ? 30.0 : 15.0;
    //             // if (segments->at (l).segment->travel_time () > 0 &&
    //             //     segments->at (l).segment->uncertainty () > 0)
    //             // {
    //             //     speed_mean = segments->at (l).segment->get_speed ();
    //             //     speed_sd = - segments->at (l).segment->length () / pow (speed_mean, 2) *
    //             //         segments->at (l).segment->uncertainty ();
    //             //     speed_sd = pow(speed_sd, 0.5);
    //             // }
    //             // else
    //             // {
    //                 // speed_mean = 10.0;
    //                 // speed_sd = 100.0;
    //             // }
    //             if (vehicle->params ()->noise_model == 0)
    //             {
    //                 double vel = speed + rng.rnorm () * vehicle->system_noise ();
    //                 while (vel <= 0 || vel > 30) vel = speed + rng.rnorm () * 3.0;//vehicle->system_noise ();
    //                 speed = vel;
    //             }
    //         }

    //         if (distance >= next_stop_d)
    //         {
    //             stop_index++;
    //             // bus arrives at stop - either stops at stop, of drives past
    //             if (bus_stop (vehicle->timestamp () - delta, rng))
    //                 distance = next_stop_d;
    //             // end of route? end.
    //             if (stop_index >= M-1)
    //             {
    //                 complete = true;
    //                 break;
    //             }
    //             // else update next stop and delta
    //             next_stop_d = stops->at (stop_index + 1).distance;
    //             delta -= (int)(dt.at (stop_index) - at.at (stop_index));
    //             if (vehicle->params ()->noise_model == 0)
    //             {
    //                 double vel = speed + rng.rnorm () * vehicle->system_noise ();
    //                 while (vel <= 0 || vel > 30) vel = speed + rng.rnorm () * 3.0;//vehicle->system_noise ();
    //                 speed = vel;
    //             }
    //             continue;
    //         }
    //     }

    //     // almost at next stop ...
    //     if (stop_index < M-1 &&
    //         next_stop_d - distance < 20)
    //     {
    //         stop_index++;
    //         bus_stop (vehicle->timestamp () - delta, rng);
    //         distance = next_stop_d;
    //         if (stop_index < M - 1)
    //             next_stop_d = stops->at (stop_index + 1).distance;
    //     }
    // }

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
        // if (u < vehicle->pr_stop ())

        /**
         * let's try Thomas' idea of stopping probability being
         * inversely proportional to speed
         * [0, 30] -> [0.9, 0.1]
         */
        double pr = (30 - speed) / 30 * 0.8 + 0.1;
        if (u < pr)
        {
            // bus stops
            // dwell = vehicle->gamma () - vehicle->dwell_time () * log (rng.runif ());
            double dwell = -1;
            while (dwell < 0)
            {
                dwell = vehicle->dwell_time () + vehicle->dwell_time_var () * rng.rnorm();
            }
            dwell += vehicle->gamma ();
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
        if (complete) return false;
        if (delta > 0) return true;

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

    void
    Particle::calculate_likelihood (latlng& y,
                                    std::vector<ShapePt>& path,
                                    double sigma)
    {
        latlng ppos = vehicle->trip ()->shape ()->coordinates_of (distance);
        double d = distanceEarth (ppos, vehicle->position ());
        log_likelihood = - log(2.0) - 2.0 * log(sigma) - pow(d, 2) / 2.0 / pow(sigma, 2);

#if VERBOSE == 2
        if (vehicle->get_n () < 20) std::cout << " => d(h(x), y) = " << d << "m";
#endif
    }

    void Particle::calculate_likelihood (Event& e, double error)
    {
        uint64_t t = (e.type == EventType::arrival) ?
            get_arrival_time (e.stop_index) :
            get_departure_time (e.stop_index);

        if (t == 0)
        {
#if VERBOSE > 0
            if (vehicle->get_n () < 20)
                std::cout << " => hasn't "
                    << (e.type == EventType::arrival ? "arrived" : "departed");
#endif
            // log_likelihood = 0.0;
            return;
        }

        int tdiff = t - e.timestamp;
#if VERBOSE > 0
        if (vehicle->get_n () < 20)
            std::cout << " => "
                << (e.type == EventType::arrival ? "arrived" : "departed")
                << " at " << t
                << " (diff: " << (tdiff) << "s)";
#endif

        log_likelihood = - 0.5 * log (2 * M_PI) - log (error) -
            pow (tdiff, 2) / 2.0 / pow (error, 2);

    }

    void Particle::calculate_arrival_likelihood (int index, uint64_t time, double error)
    {
        // particle's arrival time at stop INDEX
        if (get_arrival_time (index) == 0) log_likelihood = 0;
        log_likelihood = - 0.5 * log (2 * M_PI) - log (error) - pow (get_arrival_time (index) - time, 2) / 2.0 / pow (error, 2);
    }
    void Particle::calculate_departure_likelihood (int index, uint64_t time, double error)
    {
        if (get_departure_time (index) == 0) log_likelihood = 0;
        log_likelihood = - 0.5 * log (2 * M_PI) - log (error) - pow (get_departure_time (index) - time, 2) / 2.0 / pow (error, 2);
    }

    void Particle::set_weight (double w)
    {
        weight = w;
    }

}; // namespace Gtfs
