/**
 * Define the ETA methods.
 */

#include "etas.h"
#if VERBOSE == 2
#include "timing.h"
#endif

#ifndef NORMALAPPROX
#define NORMALAPPROX 0
#endif

#include <gsl/gsl_rng.h>
#include "vendor/rtnorm/rtnorm.hpp"

// template <typename T>
// T pnorm (T x, T m, T s)
// {
//     return 0.5 * (1.0 + erf ((x - m) * M_SQRT1_2 / s));
// }

namespace Gtfs {

    void Trip::update (uint64_t& t, RNG& rng)
    {
        // Update trip state
        if (!loaded) load ();
        bool log (false);
#if VERBOSE > 1
        log = true;
#endif

        if (log && _vehicle != nullptr)
        {
            std::ofstream fout;
            fout.open ("trip_vehicle_states.csv", std::ostream::app);
            fout << "\n" << _vehicle->vehicle_id ()
                << "," << _vehicle->trip ()->trip_id ()
                << "," << _vehicle->trip ()->route ()->route_id ()
                << "," << _vehicle->timestamp ()
                << "," << _vehicle->distance ()
                << "," << _vehicle->speed ();
            fout.close ();
        }

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
                    _segment_progress = fmax (0.0, d - Dl) / Ll;
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
                _delta = (_timestamp > 0 ? _vehicle->timestamp () - _timestamp : 0);
                _timestamp = _vehicle->timestamp ();
            }
            else
            {
                _delta = 0;
            }
        }
        else
        {
            _delta = (_timestamp > 0 ? t - _timestamp : 0);
            _timestamp = t;
            _delta = 0;
        }
        if (_delta == 0)
        {
            if (log) std::cout << " - no update, skipping.";
            return;
        }

        if (log && _vehicle != nullptr)
            std::cout
                << "\n   - Event Type: " << _event_type
                << "\n   - Stop index: " << (_stop_index+1)
                << " of " << _stops.size ()
                << "\n   - Vehicle distance: " << _vehicle->distance () << "m"
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
        if (_vehicle == nullptr) return;
        if (_delta == 0) return;

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


        gsl_rng_env_setup();                          // Read variable environnement
        const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
        gsl_rng *gen = gsl_rng_alloc (type);          // Rand generator allocation


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
                    << std::setw (4) << round (seg.segment->speed ()) << "m/s speed ("
                    << std::setw (5) << round (seg.segment->uncertainty ()) << "), "
                    << std::setw (4)
                    << (seg.segment->length () /
                        (_stops.at (i).arrival_time - _stops.at (i-1).departure_time)) << "m/s scheduled, "
                    << std::setw (4) << round (seg.segment->prior_speed ()) << " ("
                    << std::setw (5) << round (seg.segment->prior_speed_var ()) << "]";

                i++;
            }
            std::cout << "\n  ===================\n\n";
#endif


            // iterate over vehicle state
            int N (_vehicle->state ()->size ());
            N = std::min (N, 200);
            _eta_matrix = Eigen::MatrixXi::Zero (N, M);

            Particle* p;
            int l;
            int tt;
            int pi;
            double u;
            double wti;
            double dt;
            double v, xx, seg_prog;
            //  vehicle speed too slow?
            // bool use_particle_speed (true);//_vehicle->speed () < 2);
            for (int i=0; i<N; i++)
            {
                tt = 0.0;
                // fetch a random particle
                u = rng.runif ();
                wti = 0;
                pi = 0;
                while (wti < u) {
                    wti += _vehicle->state ()->at (pi++).get_weight ();
                    if (pi == _vehicle->state ()->size ()-1)
                    {
                        break;
                    }
                }
                // pi = floor (rng.runif () * (_vehicle->state ()->size ()));
                p = &(_vehicle->state ()->at (pi));
                // time to end of current segment
                l = find_segment_index (p->get_distance (), &segs);
                // std::cout << "\n --- [p], d = " << p->get_distance ()
                //     << ", l = " << l
                //     << ", D[l] = " << segs.at (l).distance;

                /**
                 * START: % of the way through segment l
                 * 0---- ... ----(l-1)-----(l)------(l+1)---- ... ----(L-1)---(L)
                 *                         [........)
                 * - if AT node (l):
                 *   - include additional dwell time
                 *   - OR wait until layover departure (whichever is later)
                 * - travel time to (l+1) -> arrival time at stop (l+1)
                 *   - if % is zero, generate a travel time from NW state
                 *   - if > zero, use particle's speed + noise OR NW state speed (if particle going slow)
                 *     or scheduled speed if NW state isn't available/useful
                 * - THEN, while l < L
                 *   - add MAX(dwell time, layover departure)
                 *   - add travel time [l, l+1) (again, NW or schedule, depending)
                 *   -> save arrival time at l+1
                 *
                 */

                /** if less than 20% of the way through the segment,
                  * or particle speed is less than 1,
                  * sample speed from segment state
                  **/
                seg_prog = (p->get_distance () - segs.at (l).distance);

                // at stop ?
                bool is_layover (_stops.at (l).departure_time > _stops.at (l).arrival_time);
                if (seg_prog < 1.0 || is_layover || l == 0)
                {
                    dt = 0.;
                    if (rng.runif () < _vehicle->pr_stop ())
                    {
                        if (_stops.at (l).sd_delay > 0)
                        {
                            dt = rtnorm (gen, 0., 5 * 60.,
                                _stops.at (l).average_delay,
                                _stops.at (l).sd_delay
                            ).first;
                        }
                        else
                        {
                            dt = rtnorm (gen, 0., 5 * 60.,
                                _vehicle->dwell_time (),
                                pow (_vehicle->dwell_time_var (), 0.5)
                            ).first;
                        }
                        dt += _vehicle->gamma ();
                    }
                    // add dwell time to travel time to NEXT stop (l+1)
                    tt += fmax (0.0, dt);

                    // first stop?
                    if (l == 0 && rng.runif () < 0.9)
                    {
                        tt = fmax (
                            0.,
                            _stops.at (l).departure_time -
                                Time (_timestamp) +
                                rng.rnorm () * 5.
                        );
                    }
                    // layover ? and Pr(driver stops) = .7
                    else if (is_layover && rng.runif () < 0.7)
                    {
                        tt = fmax (
                            tt,
                            _stops.at (l).departure_time - Time (_timestamp) + rng.rnorm () * 5.0
                        );
                    }
                }

                // std::cout << ", prog = " << seg_prog;
                double vel = 0.0;

                /** figuring out particle's speed:
                 * if there's less than 500m, OR segment progress is less than 20%,
                 * use segment state.
                 * Otherwise, sample particle's speed.
                 */

                // if (seg_prog < 100. || //(segs.at (l).segment->length () - seg_prog) > 500. ||
                //     seg_prog / segs.at (l).segment->length () < 0.1 || true)

                // more than 100m remaining? use segment state:
                if (segs.at (l).segment->length () - seg_prog > 100.)
                {
                    v = segs.at (l).segment->sample_speed (rng);
                }
                else
                {
                    auto vt = rtnorm (gen, 0.5, segs.at (l).segment->max_speed (), p->get_speed (), 2.);
                    v = vt.first;
                }
                // std::cout << ", v = " << v;
                // std::cout
                //     << " -> d = " << (segs.at (l).segment->length () - seg_prog)
                //     << " -> t = " << ((segs.at (l).segment->length () - seg_prog) / v);
                tt += (segs.at (l).segment->length () - seg_prog) / v;

                _eta_matrix (i, l+1) = tt; // no slower than 0.5m/s

                l++;
                double rho, X, Y, x, sig1, sig2, mean, var, min_tt, max_tt;
                X = 0.0;
                Y = 0.0;
                sig1 = 0.0;
                sig2 = 0.0;
                x = tt;
                while (l < L)
                {
                    dt = 0.;
                    if (rng.runif () < _vehicle->pr_stop ())
                    {
                        if (_stops.at (l).sd_delay > 0)
                        {
                            dt = rtnorm (gen, 0., 5 * 60.,
                                _stops.at (l).average_delay,
                                _stops.at (l).sd_delay
                            ).first;
                        }
                        else
                        {
                            dt = rtnorm (gen, 0., 5 * 60.,
                                _vehicle->dwell_time (),
                                pow (_vehicle->dwell_time_var (), 0.5)
                            ).first;
                        }
                        dt += _vehicle->gamma ();
                    }
                    tt += fmax (0.0, dt);

                    // layover ?
                    if (_stops.at (l).departure_time > _stops.at (l).arrival_time &&
                        rng.runif () < 0.7)
                    {
                        tt = fmax (tt, _stops.at (l).departure_time - Time (_timestamp));
                    }

                    // if (true)
                    {
                        X = Y;
                        sig1 = sig2;
                        Y = segs.at (l).segment->speed ();
                        sig2 = segs.at (l).segment->uncertainty ();
                        rho = 0.8;

                        if (sig2 > 0)
                        {
                            sig2 += pow (_delta * gtfs->parameters ()->system_noise, 2.);
                            sig2 = pow (sig2, 0.5);
                        }
                        else
                        {
                            Y = segs.at (l).segment->prior_speed ();
                            sig2 = pow (segs.at (l).segment->prior_speed_var (), 0.5);
                            if (Y == 0 || sig2 == 0)
                            {
                                Y = segs.at (l).segment->length () /
                                    (_stops.at (l+1).arrival_time - _stops.at (l).departure_time);
                                sig2 = 3.;
                            }
                        }

                        if (false) //X > 0 && sig1 > 0)
                        {
                            mean = Y + sig2 / sig1 * rho * (x - X);
                            // var(Y|X) = (1-rho^2) * sig2^2
                            var = (1 - pow (rho, 2)) * pow (sig2, 2);
                        }
                        else
                        {
                            mean = Y;
                            var = pow (sig2, 2);
                        }

                        // don't let variance get bigger than prior variance ...
                        var = fmin (var, segs.at (l).segment->prior_speed_var ());

                        // between vehicle variance
                        var += segs.at (l).segment->state_var ();

                        x = 0;
                        int n=100;
                        // x is the speed
                        auto vt = rtnorm (gen, 0.5, segs.at (l).segment->max_speed (), mean, pow (var, 0.5));
                        x = vt.first;
                        tt += segs.at (l).segment->length () / x;
                    }

                    // std::cout << "\n avg speed = ("
                    //     << segs.at (l).distance << " - " << p->get_distance ()
                    //     << ") / " << tt << " = "
                    //     << ( (segs.at (l).distance - p->get_distance ()) / tt)
                    //     << "m/s ...";

                    _eta_matrix (i, l+1) = tt;
                    l++;
                }
                // std::cout << "---\n";
            }
            // std::cout << "\n" << seg_tt.format (intMat) << "\n";

            // Convert segment mat to link mat
            // Eigen::MatrixXi link_tt = seg_tt * Hseg;

            // std::cout << "\n" << link_tt.format (intMat) << "\n";
            // std::cout << "\n" << dwell_t.format (intMat) << "\n";

            // _eta_matrix = link_tt + dwell_t;
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

        double tt_mean, tt_var, tt_lower, tt_median, tt_upper;
        double X, P, y, S, K;

#if NORMALAPPROX==0
        int normal_n;
        std::vector<double> normal_pi, normal_mu, normal_sigma2;
        std::vector<double> normal_mean, normal_var, normal_lower, normal_median, normal_upper;
        double segr; // segment proportion remaining
        // normal approx:
        // need to calculate pi, mu, sigma for each mode:
        normal_pi.clear ();
        normal_mu.clear ();
        normal_sigma2.clear ();
        normal_mean.resize (M, 0.0);
        normal_var.resize (M, 0.0);
        normal_lower.resize (M, 0.0);
        normal_median.resize (M, 0.0);
        normal_upper.resize (M, 0.0);

        if (_vehicle != nullptr && _stop_index < M - 1)
        {
            if (_stops.at (_stop_index).distance + 20.0 >= _vehicle->distance ())
            {
                // vehicle has not yet departed current stop index
                normal_pi.push_back (_vehicle->pr_stop ());
                normal_pi.push_back (_vehicle->pr_stop ());
                normal_mu.push_back (_stops.at (_stop_index).average_delay);
                normal_mu.push_back (0.0);
                normal_sigma2.push_back (pow (_stops.at (_stop_index).sd_delay, 2));
                normal_sigma2.push_back (0.0);
            }
            else
            {
                segr = _stops.at (_stop_index+1).distance - _vehicle->distance ();
                segr /= _stops.at (_stop_index+1).distance - _stops.at (_stop_index).distance;

            }
        }
        else
        {
            return etas;
        }

        // for optimisation
        const int double_bits = std::numeric_limits<double>::digits;
        std::pair<double, double> r;

        double normal_speed = 0., normal_speed_var = 0.,
            normal_tt = 0., normal_tt_var = 0.;
        auto segs = _shape->segments ();
        bool is_layover (false);
        // std::cout << "\n ++ starting from stop " << _stop_index;
        for (int m=_stop_index+1; m<M; m++)
        {
            is_layover = _stops.at (m).departure_time > _stops.at (m).arrival_time;
            // std::cout << "\n\n Estimating for stop " << m << ": ";
            // travel time to get there:
            if (segs.at (m-1).segment->uncertainty () > 0.)
            {
                normal_speed_var =
                    fmin (
                        segs.at (m-1).segment->prior_speed_var (),
                        segs.at (m-1).segment->uncertainty () +
                            normal_tt * pow(segs.at (m-1).segment->system_noise (), 2.)
                    );
                normal_speed = segs.at (m-1).segment->speed ();
            }
            else
            {
                normal_speed = segs.at (m-1).segment->prior_speed ();
                normal_speed_var = segs.at (m-1).segment->prior_speed_var ();
            }
            normal_speed_var += segs.at (m-1).segment->state_var ();

            double L = segs.at (m-1).segment->length ();
            if (m == _stop_index)
            {
                L *= segr;
            }
            normal_tt = segs.at (m-1).segment->length () / normal_speed;

            /* delta method:
             * Var(g(X)) = g'(Xhat)^2 * Var(X)
             */
            normal_tt_var = pow (segs.at (m-1).segment->length () / normal_speed, 2.);
            normal_tt_var *= normal_speed_var;
            // std::cout << " tt = " << normal_tt << " (" << normal_tt_var << ")";

            int curlen = normal_pi.size ();
            if (curlen == 0)
            {
                // initialize it
                for (int k=0;k<2;k++)
                {
                    int pk (_vehicle->pr_stop ());
                    if (k == 1) pk = 1. - pk;
                    normal_pi.push_back (pk);
                    normal_mu.push_back (_stops.at (m).average_delay * (k == 1));
                    normal_sigma2.push_back (pow (_stops.at (m).sd_delay, 2) * (k == 1));
                }
                normal_tt = 0.0;
                normal_tt_var = 0.0;
                continue;
            }
            // else if (is_layover)
            // {
            //     // driver doesn't wait, doesn't stop
            //     int pr_driver_waits = 0.7;
            //     for (int j=0;j<curlen;j++)
            //     {
            //         normal_pi.push_back (
            //             normal_pi.at (j) *
            //                 (1. - pr_driver_waits) *
            //                 (1. - _vehicle->pr_stop ())
            //         );
            //         normal_mu.push_back (normal_mu.at (j));
            //         normal_sigma2.push_back (normal_sigma2.at (j));
            //     }

            //     // driver doesn't wait, but stops
            //     for (int j=0;j<curlen;j++)
            //     {
            //         normal_pi.push_back (
            //             normal_pi.at (j) *
            //                 (1. - pr_driver_waits) *
            //                 _vehicle->pr_stop ()
            //         );
            //         normal_mu.push_back (
            //             normal_mu.at (j) + _stops.at (m).average_delay
            //         );
            //         normal_sigma2.push_back (
            //             normal_sigma2.at (j) + pow (_stops.at (m).sd_delay, 2);
            //         );
            //     }

            //     // driver waits
            //     for (int j=0;j<curlen;j++)
            //     {
            //         normal_pi.at (j) *= pr_driver_waits;
            //         int time_to_wait (_stops.at (m).departure_time - Time (_timestamp));
            //         if (time_to_wait) pr_driver_waits = 0.;
            //         if (time_to_wait > 0)
            //         {
            //             normal_mu.at (j) = time_to_wait;
            //         }

            //         // normal_sigma2.at (j) += pow (_stops.at (m).sd_delay, 2);
            //     }
            // }
            else
            {
                // it doesn't stop
                for (int j=0;j<curlen;j++)
                {

                    normal_pi.push_back (normal_pi.at (j) * (1 - _vehicle->pr_stop ()));
                    normal_mu.push_back (normal_mu.at (j));
                    normal_sigma2.push_back (normal_sigma2.at (j));
                }

                // it stops
                for (int j=0;j<curlen;j++)
                {
                    normal_pi.at (j) *= _vehicle->pr_stop ();
                    normal_mu.at (j) += _stops.at (m).average_delay;
                    normal_sigma2.at (j) += pow (_stops.at (m).sd_delay, 2);
                }

                for (int j=0;j<normal_pi.size ();j++)
                {
                    normal_mu.at (j) += normal_tt;
                    normal_sigma2.at (j) += normal_tt_var;
                }
            }

            // for (int i=0; i<std::min (20, (int)normal_pi.size ()); i++)
            // {
            //     std::cout << "\n " << i << ": "
            //         << normal_pi.at (i) << ", "
            //         << normal_mu.at (i) << ", "
            //         << normal_sigma2.at (i);
            // }

            normal_mean.at (m) = 0.0;
            double EVX = 0.0, VEX1 = 0.0, VEX2 = 0.0;
            for (int i=0; i<normal_pi.size (); i++)
            {
                normal_mean.at (m) += normal_pi.at (i) * normal_mu.at (i);
                EVX += normal_pi.at (i) * normal_sigma2.at (i);
                VEX1 += normal_pi.at (i) * pow (normal_mu.at (i), 2);
                VEX2 += normal_pi.at (i) * normal_mu.at (i);
            }
            normal_var.at (m) = EVX + VEX1 - pow (VEX2, 2);

            // std::cout << "\n arrival time -> ";

            // find quantiles:
            r = boost::math::tools::brent_find_minima (
                [&, normal_pi, normal_mu, normal_sigma2](double const& x)
                {
                    double r = 0.0;
                    for (int i=0;i<normal_pi.size ();i++)
                    {
                        r += normal_pi.at (i) *
                            R::pnorm (x, normal_mu.at (i), pow (normal_sigma2.at (i), 0.5), 1, 0);
                    }
                    return pow (r - 0.05, 2);
                },
                normal_mean.at (m) - 5 * pow (normal_var.at (m), 0.5),
                normal_mean.at (m) + 5 * pow (normal_var.at (m), 0.5),
                double_bits
            );
            normal_lower.at (m) = r.first;

            r = boost::math::tools::brent_find_minima (
                [&, normal_pi, normal_mu, normal_sigma2](double const& x)
                {
                    double r = 0.0;
                    for (int i=0;i<normal_pi.size ();i++)
                    {
                        r += normal_pi.at (i) *
                            R::pnorm (x, normal_mu.at (i), pow (normal_sigma2.at (i), 0.5), 1, 0);
                    }
                    return pow (r - 0.95, 2);
                },
                normal_mean.at (m) - 5 * pow (normal_var.at (m), 0.5),
                normal_mean.at (m) + 5 * pow (normal_var.at (m), 0.5),
                double_bits
            );
            normal_upper.at (m) = r.first;

            // std::cout << " -> ["
                // << normal_lower.at (m) << ", " << normal_upper.at (m) << "]";

            double q1, q2, qerr;
            q1 = R::qnorm (0.05, normal_mean.at (m), pow (normal_var.at (m), 0.5), 1, 0);
            q2 = R::qnorm (0.95, normal_mean.at (m), pow (normal_var.at (m), 0.5), 1, 0);
            qerr =
                pow (
                    pow (q1 - normal_lower.at (m), 2) + pow (q2 - normal_upper.at (m), 2),
                    0.5
                );
            // std::cout << "\n Quantiles of N(" << normal_mean.at (m) << ", "
            //     << normal_var.at (m) << "): [" << q1 << ", " << q2 << "] -> sqrt(sum(diff^2)) = "
            //     << qerr;


            if (qerr < 5.0 || normal_pi.size () > pow (2, 8))
            {
                // std::cout << "\n -> combining into single mode as difference is minimal!";
                normal_pi.resize (1);
                normal_pi.at (0) = 1.0;
                normal_mu.resize (1);
                normal_mu.at (0) = normal_mean.at (m);
                normal_sigma2.resize (1);
                normal_sigma2.at (0) = normal_var.at (m);
            }
        }
#endif
#if NORMALAPPROX==1
        std::vector<double> normal_mean, normal_var, normal_lower, normal_upper;
        {
            if (_vehicle == nullptr || _stop_index == M)
            {
                return etas;
            }

            normal_mean.resize (M, 0.);
            normal_var.resize (M, 0.);
            normal_lower.resize (M, 0.);
            normal_upper.resize (M, 0.);

            double tt = 0., var = 0.;
            auto segs = _shape->segments ();
            for (int m=_stop_index; m<M-1; m++)
            {
                if (m > _stop_index ||
                    _stops.at (m).distance + 20.0 >= _vehicle->distance ())
                {
                    // still at current stop: dwell time
                    tt += _stops.at (m).average_delay * _vehicle->pr_stop ();
                    var += pow (_stops.at (m).sd_delay * _vehicle->pr_stop (), 2);

                    // then segment travel time
                    tt += segs.at (m).segment->length () / segs.at (m).segment->speed ();
                    var +=
                        pow (segs.at (m).segment->length () / segs.at (m).segment->speed (), 2) *
                            segs.at (m).segment->uncertainty ();
                }
                else
                {
                    // only part of the segment remains
                    double dist_travelled (_vehicle->distance () - segs.at (m).distance);
                    tt += (segs.at (m).segment->length () - dist_travelled) /
                        segs.at (m).segment->speed ();
                    double pr (1. - dist_travelled / segs.at (m).segment->length ());
                    var += pow (pr, 2.) * segs.at (m).segment->uncertainty ();
                }


                normal_mean.at (m+1) = tt;
                normal_var.at (m+1) = var;
                normal_lower.at (m+1) = fmax (0., R::qnorm (0.025, tt, pow (var, 0.5), 1, 0));
                normal_upper.at (m+1) = R::qnorm (0.975, tt, pow (var, 0.5), 1, 0);
            }
        }
#endif

        uint64_t ts;
        if (log) std::cout << " - get state starts at stop " << (_stop_index+1) << "\n";
        for (int m=_stop_index;m<M;m++)
        {
            col_m = _eta_matrix.col (m);

            etas.at (m).stop_id = _stops.at (m).stop->stop_id ();
            etas.at (m).estimate = 0;
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

#if SIMULATION
            std::sort (
                tt_wt.begin (),
                tt_wt.end (),
                [](std::tuple<double, double>& a, std::tuple<double, double>& b)
                {
                    return std::get<0> (a) < std::get<0> (b);
                }
            );
            int i=0;
            double sum_wt = 0.0;
            // 5% quantile
            while (sum_wt < 0.05)
            {
                sum_wt += std::get<1> (tt_wt.at (i));
                i++;
            }
            tt_lower = round (std::get<0> (tt_wt.at (i-1)));
            // 50% quantile (median)
            while (sum_wt < 0.5)
            {
                sum_wt += std::get<1> (tt_wt.at (i));
                i++;
            }
            tt_median = round (std::get<0> (tt_wt.at (i-1)));
            // 90% quantile
            while (sum_wt < 0.9)
            {
                sum_wt += std::get<1> (tt_wt.at (i));
                i++;
            }
            if (i == N) i = i-1;
            tt_upper = round (std::get<0> (tt_wt.at (i)));

#endif

            /**
             * A note on indices:
             *
             * _eta_matrix is M length, INCLUDES the first stop.
             * etas/_eta_state are M length, so include first stop
             */

            // Write estimates to the file:
            // trip_id, stop_sequence, timestamp, vehicle_dist, stop_dist, pf_obs, pf_var, pf_lower, pf_upper, normal_mean, normal_var, normal_lower, normal_upper
#if SIMULATION
            fout << _trip_id << "," << _vehicle->vehicle_id () << "," << (m+1) << "," << _vehicle->timestamp ()
                << "," << _vehicle->distance () << "," << _stops.at (m).distance
                // particle predictions:
                << "," << tt_median << "," << tt_var << "," << tt_lower << "," << tt_upper
                // normal approx pred:
                << "," << normal_mean.at (m) << "," << normal_var.at (m) << ","
                << normal_lower.at (m) << "," << normal_upper.at (m);
#endif

            // update the estimate state
            if (true)
            {
                // update the state using the mean and var ...
                // std::cout << "\n     current = " << std::get<0> (_eta_state.at (m+1));
                X = (std::get<0> (_eta_state.at (m)) <= _timestamp) ? 0 :
                    std::get<0> (_eta_state.at (m)) - _timestamp;
                P = std::get<1> (_eta_state.at (m));

                if (X > 2*60*60)
                {
                    // reset X because it's going stupid!
                    // X = _stops.at (m).departure_time.asUNIX (_timestamp) - _timestamp;
                    X = _stops.at (m).average_delay;
                }
                if (P == 0)
                {
                    P = tt_var + pow (_stops.at (m).sd_delay, 2.);
                }

                if (log)
                    std::cout << "\n stop " << std::setw (2) << (m+1)
                        << " -> X = " << std::setw (8) << (round (X * 100) / 100)
                        << ", P = " << std::setw (8) << (round (P * 100) / 100);


// #if SIMULATION
//                 fout << _trip_id << "," << (m) << "," << _timestamp
//                     << "," << X << "," << P;
// #endif

                int avg_arr = _stops.at (m).departure_time.asUNIX (_timestamp) + _stops.at (m).average_delay;
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
                double eta_noise = 0.5;
                P += pow (eta_noise, 2.) * (1. + pow (_delta / (X + _delta), 2.0));
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

// #if SIMULATION
//                 fout << "," << X << "," << P << "," << tt_mean << "," << tt_var;
// #endif

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
                    y = tt_mean - X;
                    // tt_var += pow (_delta, 2.);
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
                _eta_state.at (m) = std::make_tuple (_timestamp + round (X), P);

                etas.at (m).estimate = std::get<0> (_eta_state.at (m));
                etas.at (m).quantiles.emplace_back (
                    0.025,
                    std::get<0> (_eta_state.at (m)) -
                        std::ceil (1.96 * pow (std::get<1> (_eta_state.at (m)), 0.5))
                );
                etas.at (m).quantiles.emplace_back (
                    0.5,
                    std::get<0> (_eta_state.at (m))
                );
                etas.at (m).quantiles.emplace_back (
                    0.975,
                    std::get<0> (_eta_state.at (m)) +
                        std::ceil (1.96 * pow (std::get<1> (_eta_state.at (m)), 0.5))
                );
#if SIMULATION

#endif
            }
            else
            {
                etas.at (m).estimate = ts + round (tt_mean);

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
                etas.at (m).quantiles.emplace_back (
                    0.025,
                    ts + round (std::get<0> (tt_wt.at (i-1)))
                );
                // 50% quantile (median)
                while (sum_wt < 0.5)
                {
                    sum_wt += std::get<1> (tt_wt.at (i));
                    i++;
                }
                etas.at (m).quantiles.emplace_back (
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
                etas.at (m).quantiles.emplace_back (
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
            << "              schedule         eta  delay  (prediction interval)";

        // std::string
        for (int i=0; i<stops ().size (); i++)
        {
            if (arrival_times.at (i).estimate == 0) continue;

            std::cout << "\n   + stop "
                << std::setw (2) << (i+1) << ": "
                << stops ().at (i).departure_time
                << (stops ().at (i).arrival_time < stops ().at (i).departure_time ? " *" : "  ")
                << "  "
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
