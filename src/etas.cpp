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
            Time* at;
            // needs initialization
            _eta_state.resize (_stops.size ());
            for (int m=0; m<_stops.size (); ++m)
            {
                at = &(m == 0 ? _stops.at (m).departure_time : _stops.at (m).arrival_time);
                // once prior variance is known, place here
                _eta_state.at (m) = std::make_tuple (at->asUNIX (t), 30);
            }
        }

        if (log)
        {
            std::cout
                << "\n   - Estimating ETAs using Model " << eta_model;
            std::cout << "\n   - Vehicle: "
                << (_vehicle == nullptr ? "none" : _vehicle->vehicle_id ());
        }

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

        if (log)
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
            // std::cout << "\n   == Network state ==";
            // int i=1;
            // for (auto seg : segs)
            // {
            //     std::cout << "\n ["
            //         << std::setw (5) << round (seg.distance) << "m into trip, "
            //         << std::setw (5) << round (seg.segment->length ()) << "m long, "
            //         << std::setw (3) << round (seg.segment->travel_time ()) << "s travel time ("
            //         << std::setw (4) << round (seg.segment->uncertainty ()) << "), "
            //         << std::setw (3) << (_stops.at (i).arrival_time - _stops.at (i-1).departure_time ) << "s scheduled]";
            //     i++;
            // }
            // std::cout << "\n  ===================\n\n";

            // iterate over vehicle state
            int N (_vehicle->state ()->size ());
            N = std::min (N, 10);
            Eigen::MatrixXi seg_tt (Eigen::MatrixXi::Zero (N, L));
            Eigen::MatrixXi dwell_t (Eigen::MatrixXi::Zero (N, M-1));

            Particle* p;
            int l;
            int tt;
            double dt;
            double v, xx;
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
                double rho, X, Y, x, sig1, sig2, mean, var;
                x = tt;
                while (l < L)
                {
                    // fetch tt from scheduled travel time ...
                    if (segs.at (l).segment->uncertainty () == 0 ||
                        segs.at (l).segment->uncertainty () > 2 * segs.at (l).segment->travel_time ())
                    {
                        tt += rng.rnorm () * pow (segs.at (l).segment->prior_travel_time_var (), 0.5) +
                            (int) (_stops.at (l+1).arrival_time - _stops.at (l).departure_time);
                    }
                    else
                    {
                        /**
                         * Add correlation of rho to the segments,
                         * (X,Y) ~ Normal([mu1, mu2], [[sig1,rho], [rho,sig2]) ->
                         * Y|X=x ~ Normal(mu2 + sig2/sig1 * rho * (x-mu1), (1-rho^2) * sig2^2)
                         */
                        if (segs.at (l-1).segment->uncertainty () > 0 &&
                            segs.at (l-1).segment->uncertainty () < 2 * segs.at (l-1).segment->travel_time ())
                        {
                            rho = 0.8;
                            X = segs.at (l-1).segment->travel_time ();
                            sig1 = pow (segs.at (l-1).segment->uncertainty (), 0.5);
                            Y = segs.at (l).segment->travel_time ();
                            sig2 = pow (segs.at (l).segment->uncertainty (), 0.5);
                        
                            mean = Y + sig2 / sig1 * rho * (x - X);
                            var = (1 - pow (rho, 2)) * pow (sig2, 2);
                        }
                        else
                        {
                            mean = segs.at (l).segment->travel_time ();
                            var = segs.at (l).segment->uncertainty ();
                        }
                        x = -1;
                        int n = 100;
                        while (x < segs.at (l).segment->length () / 30 | 
                               x > segs.at (l).segment->length() * 2)
                        {
                            x = rng.rnorm () * pow (var, 0.5) + mean;
                            if (n-- == 50)
                            {
                                mean = segs.at (l).segment->travel_time ();
                                var = segs.at (l).segment->uncertainty ();
                            }
                            if (n-- == 0)
                            {
                                x = segs.at (l).segment->travel_time ();
                            }
                        }
                        
                        // add between-vehicle noise (this is definitely segment independent)
                        x += rng.rnorm () * segs.at (l).segment->state_var ();
                        
                        tt += x;

                        // tt += segs.at (l).segment->sample_travel_time (rng, tt);
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
        for (int m=_stop_index;m<M-1;m++)
        {
            col_m = _eta_matrix.col (m);
            if (col_m.isZero ()) continue;

            etas.at (m+1).stop_id = _stops.at (m).stop->stop_id ();
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
                X = std::get<0> (_eta_state.at (m+1)) - _timestamp;
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

#if SIMULATION
                fout << _trip_id << "," << (m+1) << "," << _timestamp
                    << "," << X << "," << P;
#endif

                X = X;
                // 'volatility', how far to go?
                P = P + _delta;// + pow (X - tt_mean, 2);

#if SIMULATION
                fout << "," << X << "," << P << "," << tt_mean << "," << tt_var;
#endif

                // KF update
                if (tt_var == 0)
                {
                    tt_var = tt_mean;
                }
                else
                {
                    tt_var += tt_mean;
                    // tt_var += pow (X - tt_mean, 2);
                }
                y = tt_mean - X;
                S = P + tt_var;
                K = P / S;
                X += K * y;
                P = (1 - K) * P;

                P += pow (y, 2);

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
                << (arrival_times.at (i).estimate - stops ().at (i).arrival_time)
                << "  ("
                << Time (arrival_times.at (i).quantiles.at (0).time)
                << " - "
                << Time (arrival_times.at (i).quantiles.at (2).time)
                << ")";
        }
    }

} // namespace Gtfs
