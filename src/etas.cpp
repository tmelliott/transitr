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

        bool log (false);
#if VERBOSE > 1
        log = true;
#endif

        if (log)
            std::cout << "\n * Updating trip " << _trip_id
                << " (" << _route->route_short_name () << ")";

        int eta_model (gtfs->parameters ()->eta_model);

        // Initalize?
        if (_eta_state.size () == 1)
        {

            /** Check if `stop_delays` table exists, and if so for each stop fetch the
             * relevant delay information.
             */

            get_db_delays ();

            Time* at;
            // needs initialization
            _eta_state.resize (_stops.size ());
            double d = 0, e = 30;
            for (int m=0; m<_stops.size (); ++m)
            {
                at = &(m == 0 ? _stops.at (m).departure_time : _stops.at (m).arrival_time);

                d = _stops.at (m).average_delay;
                e = _stops.at (m).sd_delay;

                if (e == 0) e = 600;

                // once prior variance is known, place here
                _eta_state.at (m) = std::make_tuple (at->asUNIX (t) + d, pow (e, 2));
            }
        }

        if (log)
        {
            std::cout
                << "\n   - Estimating ETAs using Model " << eta_model;
            std::cout << "\n   - Vehicle: "
                << (_vehicle == nullptr ? "none" : _vehicle->vehicle_id ());
        }

        _segment_progress = -1;
        if (_vehicle != nullptr)
        {
            // has vehicle been updated?
            if (_vehicle->timestamp () > _timestamp)
            {
                if (log) std::cout << " (updated - ";

                Event* event = _vehicle->latest_event ();

                _segment_progress = 0;
                if (event->type == EventType::gps)
                {
                    if (log) std::cout << "gps";
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
                        if (log) std::cout << "arrival";
                    }
                    else
                    {
                        _event_type = 2;
                        if (log) std::cout << "departure";
                    }

                    _stop_index = event->stop_index;
                    double d = _stops.at (_stop_index).distance;
                    _segment_index = find_segment_index (d, &(_shape->segments ()));
                }
                if (log) std::cout << ")";

            }
        }
        _delta = (_timestamp > 0 ? t - _timestamp : 0);
        _timestamp = t;

        if (log && _vehicle != nullptr)
            std::cout
                << "\n   - Event Type: " << _event_type
                << "\n   - Stop index: " << (_stop_index+1)
                << " of " << _stops.size ()
                << "\n   - Segment index: " << (_segment_index+1)
                << " of " << _shape->segments ().size ()
                << " (" << (_segment_progress*100) << "%)";

        if (_stop_index < _shape->segments ().size ())
        {
            if (log) std::cout << "\n\n * Forecasting ..." << std::endl;
            forecast (rng);

            if (log) std::cout << "\n\n * Summarising ..." << std::endl;
            arrival_times = get_etas ();

            if (log) print_etas ();
        }

        if (log) std::cout << "\n\n------------------------------------------\n";

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
            // current NW state
#if VERBOSE > 2
            std::cout << "\n   == Network state ==";
            int i=1;
            for (auto seg : segs)
            {
                std::cout << "\n ["
                    << std::setw (5) << round (seg.distance) << "m into trip, "
                    << std::setw (5) << round (seg.segment->length ()) << "m long, "
                    << std::setw (4) << round (seg.segment->travel_time ()) << "s travel time ("
                    << std::setw (5) << round (seg.segment->uncertainty ()) << "), "
                    << std::setw (4) << (_stops.at (i).arrival_time - _stops.at (i-1).departure_time ) << "s scheduled, "
                    << std::setw (4) << round (seg.segment->prior_travel_time ()) << " ("
                    << std::setw (5) << round (seg.segment->prior_travel_time_var ()) << "]";

                i++;
            }
            std::cout << "\n  ===================\n\n";
#endif

            // iterate over vehicle state
            int N (_vehicle->state ()->size ());
            N = std::min (N, 500);
            Eigen::MatrixXi seg_tt (Eigen::MatrixXi::Zero (N, L));
            Eigen::MatrixXi dwell_t (Eigen::MatrixXi::Zero (N, M-1));

            Particle* p;
            int l;
            int tt;
            int pi;
            double dt;
            double v, xx, seg_prog;
            //  vehicle speed too slow?
            // bool use_particle_speed (true);//_vehicle->speed () < 2);
            for (int i=0; i<N; i++)
            {
                tt = 0;
                // fetch a random particle
                pi = floor (rng.runif () * (_vehicle->state ()->size ()));
                p = &(_vehicle->state ()->at (pi));
                // time to end of current segment
                l = find_segment_index (p->get_distance (), &segs);

                /** if less than 20% of the way through the segment,
                  * or particle speed is less than 1,
                  * sample speed from segment state
                  **/
                seg_prog = (p->get_distance () - segs.at (l).distance) / segs.at (l).segment->length ();
                if (_vehicle->speed () > 2.0 &&
                    seg_prog > 0.2 &&
                    p->get_speed () > 1.0 &&
                    rng.runif () < 0.5)
                {
                    v = fmin (30.0, fmax (1.0, rng.rnorm () * 3.0 + p->get_speed ()));
                }
                else
                {
                    v = segs.at (l).segment->sample_speed (rng);
                }

                xx = segs.at (l).segment->sample_travel_time (rng);
                tt += fmin(xx,
                    round ((
                        segs.at (l).distance +
                        segs.at (l).segment->length () -
                        p->get_distance ()
                    ) / p->get_speed ())
                );

                seg_tt (i, l) = tt; // no slower than 0.5m/s

                l++;
                double rho, X, Y, x, sig1, sig2, mean, var, min_tt, max_tt;
                X = 0.0;
                Y = 0.0;
                sig1 = 0.0;
                sig2 = 0.0;
                x = tt;
                while (l < L)
                {
                    // Q: is the next stop a layover?
                    if (_stops.at (l).departure_time > _stops.at (l).arrival_time)
                    {
                        // time until scheduled departure
                        tt = fmax (_stops.at (l).departure_time - Time (_timestamp), tt);
                    }

                    // fetch tt from scheduled travel time ...
                    // if (segs.at (l).segment->uncertainty () == 0 ||
                    //     segs.at (l).segment->uncertainty () > 2 * segs.at (l).segment->travel_time ())
                    // {
                    //     if (segs.at (l).segment->prior_travel_time_var () == 0)
                    //     {
                    //         // travel time is just going to be somewhere between 5m/s and 30m/s
                    //         double len = segs.at (l).segment->length ();
                    //         double sched_speed = len /
                    //             ((int) (_stops.at (l+1).arrival_time - _stops.at (l).departure_time));
                    //         double v = 0;
                    //         int n = 100;
                    //         while (v < 2 || v > 30)
                    //         {
                    //             v = rng.rnorm () * 5.0 + sched_speed;

                    //             // give up after 100 tries
                    //             if (n-- == 0)
                    //                 v = rng.runif () * 28.0 + 2.0;
                    //         }
                    //         tt += len / v;
                    //     }
                    //     else
                    //     {
                    //         // prior uncertainty + vehicle variation
                    //         double var = segs.at (l).segment->prior_travel_time_var ();
                    //         // var += pow (segs.at (l).segment->state_var (), 2);
                    //         double tx = 0;
                    //         double len = segs.at (l).segment->length ();
                    //         int nn = 100;
                    //         while (tx < len / 30 || tx > len)
                    //         {
                    //             tx = rng.rnorm () * pow (var, 0.5) +
                    //                 (int) (_stops.at (l+1).arrival_time - _stops.at (l).departure_time);
                    //             if (nn-- == 0)
                    //             {
                    //                 var = pow (len, 0.5);
                    //             }
                    //         }
                    //         tt += tx;
                    //     }
                    // }
                    // else
                    // {
                    //     /**
                    //      * Add correlation of rho to the segments,
                    //      * (X,Y) ~ Normal([mu1, mu2], [[sig1,rho], [rho,sig2]) ->
                    //      * Y|X=x ~ Normal(mu2 + sig2/sig1 * rho * (x-mu1), (1-rho^2) * sig2^2)
                    //      */
                    //     if (segs.at (l-1).segment->uncertainty () > 0 &&
                    //         segs.at (l-1).segment->uncertainty () < 2 * segs.at (l-1).segment->travel_time ())
                    //     {
                    //         rho = 0.0;
                    //         X = segs.at (l-1).segment->travel_time ();
                    //         sig1 = pow (segs.at (l-1).segment->uncertainty (), 0.5);
                    //         Y = segs.at (l).segment->travel_time ();
                    //         sig2 = pow (segs.at (l).segment->uncertainty (), 0.5);

                    //         mean = Y + sig2 / sig1 * rho * (x - X);
                    //         var = (1 - pow (rho, 2)) * pow (sig2, 2);
                    //     }
                    //     else
                    //     {
                    //         mean = segs.at (l).segment->travel_time ();
                    //         var = segs.at (l).segment->uncertainty ();
                    //     }
                    //     x = -1;
                    //     int n = 100;
                    //     while (x < segs.at (l).segment->length () / 30 |
                    //            x > segs.at (l).segment->length() * 2)
                    //     {
                    //         x = rng.rnorm () * pow (var, 0.5) + mean;
                    //         if (n-- == 50)
                    //         {
                    //             mean = segs.at (l).segment->travel_time ();
                    //             var = segs.at (l).segment->uncertainty ();
                    //         }
                    //         if (n-- == 0)
                    //         {
                    //             x = segs.at (l).segment->travel_time ();
                    //         }
                    //     }

                    //     // add between-vehicle noise (this is definitely segment independent)
                    //     double newx;
                    //     while (newx < segs.at (l).segment->length () / 30 |
                    //            newx > segs.at (l).segment->length() * 2)
                    //     {
                    //         newx = x + rng.rnorm () * segs.at (l).segment->state_var ();
                    //     }
                    //     tt += newx;

                    //     // tt += segs.at (l).segment->sample_travel_time (rng, tt);
                    // }

                    {
                        // X = Y;
                        // sig1 = sig2;
                        Y = segs.at (l).segment->travel_time ();
                        sig2 = pow (segs.at (l).segment->uncertainty (), 0.5);

                        rho = 0.8;

                        if (sig2 == 0 || sig2 > 2 * Y)
                        {
                            Y = segs.at (l).segment->prior_travel_time ();
                            sig2 = pow (segs.at (l).segment->prior_travel_time_var (), 0.5);
                            if (Y == 0 || sig2 == 0)
                            {
                                Y = _stops.at (l+1).arrival_time - _stops.at (l).departure_time;
                                sig2 = pow (3.0 + 0.3 * Y, 0.5);
                            }
                        }

                        if (X > 0 && sig1 > 0)
                        {
                            mean = Y + sig2 / sig1 * rho * (x - X);
                            var = (1 - pow (rho, 2)) * pow (sig2, 2);
                        }
                        else
                        {
                            mean = Y;
                            var = pow (sig2, 2);
                        }


                        min_tt = segs.at (l).segment->length () / 30.0;
                        max_tt = segs.at (l).segment->length () * 2.0;

                        // between vehicle variance
                        var += segs.at (l).segment->state_var ();

                        x = 0;
                        int n=100;
                        while (x < min_tt || x > max_tt)
                        {
                            x = rng.rnorm () * pow (var, 0.5) + mean;
                            if (n-- == 0)
                            {
                                x = rng.runif () * (max_tt/2 - min_tt) + min_tt;
                            }
                        }

                        tt += x;
                    }
                    seg_tt (i, l) = tt;
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
            // std::cout << "\n" << seg_tt.format (intMat) << "\n";

            // Convert segment mat to link mat
            Eigen::MatrixXi link_tt = seg_tt * Hseg;

            // std::cout << "\n" << link_tt.format (intMat) << "\n";
            // std::cout << "\n" << dwell_t.format (intMat) << "\n";

            _eta_matrix = link_tt + dwell_t;
            // std::cout << "\n" << _eta_matrix.format (intMat) << "\n";
        }
        else
        {

        }
    }

    etavector Trip::get_etas ()
    {
        bool log (false);
#if VERBOSE > 1
        log = true;
#endif
        etavector etas;
        int M (_stops.size ());
        etas.resize (M);
        int N (_eta_matrix.rows ());

        Eigen::VectorXi col_m;
        std::vector<std::tuple<double, double> > tt_wt (N);

#if SIMULATION
        // Want to write them to a file (for research)
        std::ofstream fout;
        fout.open ("eta_state.csv", std::ofstream::app);
        // trip_id, stop_id, timestamp, state, var, pred_state, pred_var, obs, obs_err, final_state, final_var
#endif

        double tt_mean, tt_var, tt_lower, tt_upper;
        double X, P, y, S, K;
        uint64_t ts;
        // std::cout << " - get state starts at stop " << _stop_index << "\n";
        for (int m=_stop_index;m<M-1;m++)
        {
            col_m = _eta_matrix.col (m);

            etas.at (m+1).stop_id = _stops.at (m).stop->stop_id ();
            etas.at (m+1).estimate = 0;
            if (col_m.isZero ()) continue;

            for (int i=0; i<N; ++i)
            {
                if (_vehicle != nullptr && _vehicle->state()->size () == N)
                {
                    tt_wt.at (i) = std::make_tuple (
                        col_m (i),
                        _vehicle->state ()->at (i).get_weight ()
                    );
                }
                else
                {
                    tt_wt.at (i) = std::make_tuple (col_m (i), pow(N, -1));
                }
            }

            tt_mean = std::accumulate (tt_wt.begin (), tt_wt.end (), 0.0,
                [](double a, std::tuple<double, double>& b)
                {
                    if (std::get<0> (b) == 0) return a;
                    return a + std::get<0> (b) * std::get<1> (b);
                });

            tt_var = std::accumulate (tt_wt.begin (), tt_wt.end (), 0.0,
                [&tt_mean](double a, std::tuple<double, double>& b)
                {
                    return a + pow (std::get<0> (b) - tt_mean, 2) * std::get<1> (b);
                });

            ts = (_vehicle != nullptr && _vehicle->state ()->size () == N) ?
                _vehicle->timestamp () : _timestamp;

            /**
             * A note on indices:
             *
             * _eta_matrix is M-1 length, so EXCLUDES the first stop.
             * etas/_eta_state are M length, so include first stop
             */

            // update the estimate state
            if (true)
            {
                // update the state using the mean and var ...
                // std::cout << "\n     current = " << std::get<0> (_eta_state.at (m+1));
                X = (std::get<0> (_eta_state.at (m+1)) <= _timestamp) ? 0 :
                    std::get<0> (_eta_state.at (m+1)) - _timestamp;
                P = std::get<1> (_eta_state.at (m+1));

                if (X > 2*60*60)
                {
                    // reset X because it's going stupid!
                    X = _stops.at (m+1).departure_time.asUNIX (_timestamp) - _timestamp;
                }
                if (P == 0)
                {
                    P = abs (X);
                }

                if (log)
                    std::cout << "\n stop " << std::setw (2) << m
                        << " -> X = " << std::setw (8) << (round (X * 100) / 100)
                        << ", P = " << std::setw (8) << (round (P * 100) / 100);


#if SIMULATION
                fout << _trip_id << "," << (m+1) << "," << _timestamp
                    << "," << X << "," << P;
#endif

                int avg_arr = _stops.at (m+1).departure_time.asUNIX (_timestamp) + _stops.at (m+1).average_delay;
                double avg_eta = avg_arr - _timestamp;
                // if (_stops.at (m+1).sd_delay > 0.5 && avg_eta > 0.0)
                // {
                //     double F, lambda;
                //     lambda = 0.5; // ratio of prior variance and how long to go?
                //     F = lambda * (X - avg_eta) / X;
                //     std::cout << "\n F = " << F << ", X = " << X << ", avg = " << avg_eta;
                //     X = F * X;
                //     P = F * P * F;
                // }
                P += _delta;
                // if (avg_eta > 0) P += avg_eta;

                // X = F * X;
                // // 'volatility', how far to go?
                // P = F * P * F + _delta;// + pow (X - tt_mean, 2);

                // tt_var = fmax (pow (tt_mean, 2), tt_var);
                // tt_var += tt_mean * 2;

                // multiple by how far it is from state estimate
                // tt_var *= fmax (1.0, (tt_mean - X) / pow (P, 0.5));

                // if (tt_var < 2 * X)
                // {
                //     tt_var = pow (20 * sqrt (fmin (10, tt_mean / 60)) + 30, 2);
                // }
                // else
                // {
                //     tt_var = 1e10; // very, very unreliable!
                // }

                if (log)
                    std::cout
                        << " => Z = " << std::setw (8) << (round (tt_mean * 100) / 100)
                        << ", R = " << std::setw (8) << (round (tt_var * 100) / 100);

#if SIMULATION
                fout << "," << X << "," << P << "," << tt_mean << "," << tt_var;
#endif

                // KF update
                // if (tt_var == 0)
                // {
                //     tt_var = tt_mean;
                // }
                // else
                // {
                //     tt_var += tt_mean;
                //     // tt_var += pow (X - tt_mean, 2);
                // }

                if (tt_var < 1e6 && tt_mean < 2*60*60)
                {
                    tt_var += pow (tt_mean, 2);
                    y = tt_mean - X;
                    S = P + tt_var;
                    K = P / S;
                    X += K * y;
                    P = (1 - K) * P;
                }

                if (log)
                    std::cout
                        << " -> X = " << std::setw (8) << (round (X * 100) / 100)
                        << ", P = " << std::setw (8) << (round (P * 100) / 100);

                // P += pow (y, 2);

#if SIMULATION
                fout << "," << X << "," << P << "\n";
#endif

                // insert back into state
                _eta_state.at (m+1) = std::make_tuple (_timestamp + round (X), P);

                etas.at (m+1).estimate = std::get<0> (_eta_state.at (m+1));
                etas.at (m+1).quantiles.emplace_back (
                    0.025,
                    std::get<0> (_eta_state.at (m+1)) -
                        std::ceil (1.96 * pow (std::get<1> (_eta_state.at (m+1)), 0.5))
                );
                etas.at (m+1).quantiles.emplace_back (
                    0.5,
                    std::get<0> (_eta_state.at (m+1))
                );
                etas.at (m+1).quantiles.emplace_back (
                    0.975,
                    std::get<0> (_eta_state.at (m+1)) +
                        std::ceil (1.96 * pow (std::get<1> (_eta_state.at (m+1)), 0.5))
                );
            }
            else
            {
                etas.at (m+1).estimate = ts + round (tt_mean);

                std::sort (tt_wt.begin (), tt_wt.end (),
                    [](std::tuple<double, double>& a, std::tuple<double, double>& b)
                    {
                        return std::get<0> (a) < std::get<0> (b);
                    });

                // then the quantiles
                int i=0;
                double sum_wt = 0.0;
                // 2.5% quantile
                while (sum_wt < 0.025)
                {
                    sum_wt += std::get<1> (tt_wt.at (i));
                    i++;
                }
                etas.at (m+1).quantiles.emplace_back (
                    0.025,
                    ts + round (std::get<0> (tt_wt.at (i-1)))
                );
                // 50% quantile (median)
                while (sum_wt < 0.5)
                {
                    sum_wt += std::get<1> (tt_wt.at (i));
                    i++;
                }
                etas.at (m+1).quantiles.emplace_back (
                    0.5,
                    ts + round (std::get<0> (tt_wt.at (i-1)))
                );
                // 97.5% quantile
                while (sum_wt < 0.975)
                {
                    sum_wt += std::get<1> (tt_wt.at (i));
                    i++;
                }
                if (i == N) i = i-1;
                etas.at (m+1).quantiles.emplace_back (
                    0.975,
                    ts + round (std::get<0> (tt_wt.at (i)))
                );
            }
        }

#if SIMULATION
        {
            fout.close ();
        }
#endif

        return etas;
    }

    etavector& Trip::get_arrival_times ()
    {
        return arrival_times;
    }

    void Trip::print_etas ()
    {
        // if (!state_initialised) return;

        std::cout << "\n\n - ETAs for route "
            << route ()->route_short_name ()
            << " ("
            << stops ().at (0).arrival_time
            << ") with " << stops ().size () << " stops\n"
            << "              schedule       eta  delay  (prediction interval)";

        // std::string
        for (int i=0; i<stops ().size (); i++)
        {
            if (arrival_times.at (i).estimate == 0) continue;

            std::cout << "\n   + stop "
                << std::setw (2) << (i+1) << ": "
                << stops ().at (i).arrival_time << "  "
                << Time (arrival_times.at (i).estimate) << "  "
                << std::setw (5)
                << (arrival_times.at (i).estimate - stops ().at (i).arrival_time);
            if (arrival_times.at (i).quantiles.size () == 3 &&
                arrival_times.at (i).quantiles.at (0).time > 0 &&
                arrival_times.at (i).quantiles.at (2).time > 0)
            {
                std::cout
                    << "  ("
                    << Time (arrival_times.at (i).quantiles.at (0).time)
                    << " - "
                    << Time (arrival_times.at (i).quantiles.at (2).time)
                    << ")";
            }
            else
            {
                std::cout << "  (NA)";
            }
        }
    }

} // namespace Gtfs
