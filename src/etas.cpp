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
        if (!valid () || complete () || _delta == 0) return;

        /**
         * Predict arrival times for a vehicle using network state
         * and current position.
         *
         * - network state gives full-segment travel times, which
         *   can be summed to get arrival time
         * - current position gives travel time of current segment,
         *   which is estimated from particle fleet
         */

#if VERBOSE == 2
        Timer timer;
        std::cout << "\n- vehicle " << _vehicle_id << " - predicting etas";
#endif

        auto segments = _trip->shape ()->segments ();
        auto stops = _trip->stops ();
        int M = stops.size ();
        int L = segments.size ();

        Eigen::IOFormat decimalMat (6, 0, ", ", "\n", "  [", "]");
        Eigen::IOFormat ColVec (6, 0, ", ", "\n", "  [", "]");
        Eigen::IOFormat tColVec (6, 0, " ", ", ", "", "", "  [", "]^T");

        // Convert from state to proper ETAs
        Eigen::MatrixXd H_eta (Eigen::MatrixXd::Zero (M, M));
        for (int i=0; i<M; i++) for (int j=0; j<=i; j++) H_eta (i, j) = 1;

        // Determine Z and R
        Eigen::VectorXd Z_seg (L);
        Eigen::MatrixXd R_seg (Eigen::MatrixXd::Zero (L, L));

        for (int i=0; i<L; i++)
        {
            Z_seg (i) = segments.at (i).segment->travel_time ();
            R_seg (i, i) = segments.at (i).segment->uncertainty () == 0 ?
                pow (segments.at (i).segment->travel_time (), 0.5) :
                segments.at (i).segment->uncertainty ();
        }
#if VERBOSE == 2
        std::cout << "\n\n * segment travel time estimates:";
        std::cout << "\n\n  >> Z_seg:\n" << Z_seg.format (tColVec);
        std::cout << "\n\n  >> R_seg:\n" << R_seg.format (decimalMat);
#endif

        // Convert them to STOP travel times
        Eigen::MatrixXd H_seg (Eigen::MatrixXd::Zero (M, L));
        int l = 0;
        double d0, d1;
        for (int j=1; j<M; j++)
        {
            // std::cout << "\n stop " << (j-1) << " to " << j << ": ";
            // stop to end of current segment (should be ~0)
            d0 = stops.at (j-1).distance;
            d1 = stops.at (j).distance;
            double seglen;
            while (l < L)
            {
                // length of this segment * proportion in stop
                seglen = segments.at (l).segment->length ();
                if (segments.at (l).distance < d0)
                {
                    // remove length from segment start to stop
                    seglen -= (d0 - segments.at (l).distance);
                }
                if (l < L && segments.at (l).distance + segments.at (l).segment->length () > d1)
                {
                    // remove end of segment
                    seglen -= (segments.at (l).distance + segments.at (l).segment->length () - d1);
                }
                // std::cout << " -> " << (seglen / segments.at (l).segment->length ()) << " * seg[" << l << "]";

                H_seg (j, l) = seglen / segments.at (l).segment->length ();

                if (segments.at (l).distance + segments.at (l).segment->length () > d1) break;
                l++;
            }
        }

        Eigen::VectorXd Z (M);
        Eigen::MatrixXd R (Eigen::MatrixXd::Zero (M, M));
        Z = H_seg * Z_seg;
        R = H_seg * R_seg * H_seg.transpose ();
        R (0, 0) = 1e6;
#if VERBOSE == 2
        std::cout << "\n\n * inter-stop travel time estimates:";
        std::cout << "\n\n  >> H_seg:\n" << H_seg.format (decimalMat);
        std::cout << "\n\n  >> Z:\n" << Z.format (tColVec);
        std::cout << "\n\n  >> R:\n" << R.format (decimalMat);
        std::cout << "\n\n  >> delta: " << _delta << " seconds";
#endif

        // Current ETA state:
        Eigen::VectorXd ETA_hat (M);
        Eigen::MatrixXd ETA_var (Eigen::MatrixXd::Zero (M, M));
        for (int i=0; i<M; i++) 
        {
            ETA_hat (i) = _tt_state.at (i);
            for (int j=0; j<M; j++)
            {
                ETA_var (i ,j) = _tt_cov.at (i).at (j);
            }
        }
        // ETA_hat (M) = Time (_timestamp - _delta) - stops.at (0).departure_time;
        // also need a condiiton for "bus hasn't started yet" to add data on ETA[0]

#if VERBOSE == 2
        std::cout << "\n\n * initial state";
        std::cout << "\n\n  >> E(ETA):\n" << ETA_hat.format (tColVec);
        std::cout << "\n\n  >> Var(ETA):\n" << ETA_var.format (decimalMat);  
        Eigen::VectorXd ETA_hat_original (ETA_hat);
        Eigen::MatrixXd ETA_original (ETA_var);
        std::cout << "\n\n  >> ETAs:\n" << (H_eta * ETA_hat).format (tColVec);
#endif

        // Compute F, H, Q
        Eigen::MatrixXd I (Eigen::MatrixXd::Identity (M, M));
        Eigen::MatrixXd F (I);
        Eigen::MatrixXd H (I);
        Eigen::MatrixXd Q (I);
        Q = Q * _delta;
#if VERBOSE == 2
        std::cout << "\n\n * control matrices";
        std::cout << "\n\n  >> F:\n" << F.format (decimalMat);
        std::cout << "\n\n  >> H:\n" << H.format (decimalMat);
        std::cout << "\n\n  >> Q:\n" << Q.format (decimalMat);
#endif

        // predicted state
        ETA_hat = F * ETA_hat;
        ETA_var = F * ETA_var * F.transpose () + Q;
#if VERBOSE == 2
        std::cout << "\n\n * predicted state";
        std::cout << "\n\n  >> E(ETA):\n" << ETA_hat.format (tColVec);
        std::cout << "\n\n  >> Var(ETA):\n" << ETA_var.format (decimalMat);
        std::cout << "\n\n  >> ETAs:\n" << (H_eta * ETA_hat).format (tColVec);
#endif

        // updated state
        auto y = Z - H * ETA_hat;
        auto S = R + H * ETA_var * H.transpose ();
        auto K = ETA_var * H.transpose () * S.inverse ();

        ETA_hat += K * y;
        ETA_var = (I - K * H) * ETA_var * (I - K * H).transpose () + K * R * K.transpose ();
#if VERBOSE == 2
        std::cout << "\n\n * updated state";
        std::cout << "\n\n  >> E(ETA):\n" << ETA_hat.format (tColVec);
        std::cout << "\n\n  >> Var(ETA):\n" << ETA_var.format (decimalMat);
        std::cout << "\n\n >> ETAs:\n" << (H_eta * ETA_hat).format (tColVec);
#endif

        // Now figure out how far into the trip the vehicle currently is,
        // _ignorant of previous stops which could be wrong_.
        _current_stop = find_stop_index (distance (), &stops);
#if VERBOSE == 2
        std::cout << "\n\n - current segment: " << (_current_stop + 1) << " of " << stops.size ();
#endif
        if (_current_stop == stops.size () - 1)
        {
#if VERBOSE == 2
            std::cout << " ... well that's awkward, this vehicle has finished!";
#endif
        }
        else
        {
            float pr = (distance () - stops.at (_current_stop).distance) /
                (stops.at (_current_stop + 1).distance - stops.at (_current_stop).distance);
#if VERBOSE == 2
            std::cout << "\n - proportion of segment complete: ";
            std::cout << pr;
#endif
            // now adjust ETAs ... 
            int t_passed (Time (_timestamp) - trip_start_time ());
#if VERBOSE == 2
            std::cout << "\n - time since scheduled trip start: " << t_passed << "s\n";
#endif
            if (_current_stop > 0 || t_passed > 0)
            {
                for (int i=0; i<_current_stop; i++) t_passed -= ETA_hat (i);
                t_passed -= pr * ETA_hat (_current_stop);
                ETA_hat (0) = t_passed;
            }
#if VERBOSE == 2
            std::cout << "\n\n * updated state (to get valid arrival time)";
            std::cout << "\n\n  >> E(ETA):\n" << ETA_hat.format (tColVec);
            std::cout << "\n\n  >> Var(ETA):\n" << ETA_var.format (decimalMat);
            std::cout << "\n\n >> ETAs:\n" << (H_eta * ETA_hat).format (tColVec);
#endif

        }

        // update vehicle's state
        for (int i=0; i<M; i++)
        {
            _tt_state.at (i) = ETA_hat (i);
            for (int j=0; j<M; j++)
            {
                _tt_cov.at (i).at (j) = ETA_var (i, j);
            }
        }
        _tt_time = _timestamp;

#if VERBOSE == 2
        std::cout << "\n\n ==> FINAL ETAs:";
        int tarr = this->trip_start_time ().seconds ();
        for (int i=0; i<M; i++)
        {
            tarr += ETA_hat (i);
            std::cout << "\n   - stop " << std::setw(2) << (i+1) << ": "
                << Time (tarr);
        }
#endif

        return;


#if VERBOSE == 2
        if (_N < 20) std::cout << "\n    --- Particle ETAs ...";
#endif
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

#if VERBOSE == 2
            if (_N < 20)
            {
                std::cout << "\n    => ";
                for (auto eta : p->get_arrival_times ()) {
                    if (eta <= _timestamp) std::cout << "*, ";
                    else std::cout << ((eta - _timestamp)) << ", ";
                }
            }
#endif
        }

#if VERBOSE == 2
        if (_N < 20)
        {
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
        }
#endif

        /**
         * Now, we assume the particles have taken into account any correlation
         * structure between segment travel times.
         *
         * Simply need to compute B and Cov matrix for travel times (for this vehicle)
         */
        
        // auto stops = _trip->stops ();
        // int M = stops.size ();
        
#if VERBOSE == 2
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
#endif
        std::vector<double> tt_obs (M, 0.0);
        std::vector<std::vector<double> > tt_r (M, std::vector<double> (M, 0.0));

        double cov_lj;
        // std::cout << "  > curtime = " << _timestamp << "\n";
        // std::cout << "  > curstop = " << _current_stop << "\n";
        for (int l=0; l<M; l++)
        {
#if VERBOSE == 2
            if (_N < 20) std::cout << "\n STOP " << l << ": ";
#endif
            tt_obs.at (l) = std::accumulate(_state.begin (), _state.end (), 0.0,
                                            [&](double d, Particle& p) {
                                                int dt;
                                                // std::cout << "[" << p.get_stop_index () << "]";
                                                if (l <= p.get_stop_index ())
                                                {
                                                    // this stop is behind this particle
                                                    dt = 0;
                                                }
                                                else if (l == p.get_stop_index () + 1)
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
#if VERBOSE == 2
                                                if (_N < 20) std::cout << std::setw (5) << dt << ", ";
#endif
                                                return d + dt;
                                            });
            tt_obs.at (l) /= (double)_state.size ();
        }
#if VERBOSE == 2
        std::cout << "\n\n Z vector: [";
        for (auto z : tt_obs) std::cout << " " << z << " ";
#endif

        for (int l=0; l<M; l++) 
        {
            for (int j=0; j<M; j++)
            {
                // if (l != j)
                // {
                //     tt_r.at (l).at (j) = 0;
                //     continue;
                // }
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

#if VERBOSE == 2
        std::cout << "]\n R matrix: ";
        for (auto br : tt_r)
        {
            std::cout << "\n  ";
            for (auto bc : br) std::cout << std::setw(9) << std::round (bc) << "  ";
        }
        std::cout << "\n ---------\n";

        std::cout << "\n Now update the travel time estimates' state ...\n";
#endif
        int tt_delta;
        if (_tt_time == 0)
        {
            // this here will need to be fixed up though
            double tsum = 0, rsum;
            for (int i=0; i<M; i++)
            {
                if (tt_obs.at (i) == 0)
                {

                } 
                else if (tsum == 0)
                {
                    // partial segment, so use particle speeds to set "prior"
                    tsum += tt_obs.at (i);
                }
                else
                {
                    // use scheduled time difference between stops
                    tsum += _trip->stops ().at (i).arrival_time - 
                        _trip->stops ().at (i-1).arrival_time;
                }
                _tt_state.at (i) = tsum;

                for (int j=0; j<M; j++)
                {
                    // 5-minute schedule deviation at each stop (variance)
                    _tt_cov.at (i).at (j) = (int)(i == j) * 100.0 * i;
                }
            }
            // _tt_cov = tt_r; // that can't have been good ... D=
            tt_delta = 0;
        }
        else
        {
            tt_delta = _timestamp - _tt_time;
        }

        if (tt_delta > 0)
        {
#if VERBOSE == 2
            std::cout << " -> ETA delta: " << tt_delta << " seconds ... \n";
#endif

            // we need to work on the submatrix (exclude passed stops BASED ON THE DATA (not the state))
            int min_i = 0;
            while (tt_obs.at (min_i) == 0 && min_i < M-1) min_i++;
            M -= min_i;

            // std::cout << "From stop " << min_i << " (M = " << M << ")\n";

            Eigen::IOFormat decimalMatFmt (2, 0, ", ", "\n", "[", "]");

            // --- Predict equations
            // Convert B to (column) vector
            Eigen::VectorXd Bmat (M);
            for (int i=0; i<M; i++) Bmat (i) = _tt_state.at (i+min_i);
#if VERBOSE == 2
            std::cout << "\n Bhat = \n" << Bmat;
#endif

            // Convert P to a matrix
            Eigen::MatrixXd Pmat (M, M);
            for (int i=0; i<M; i++)
                for (int j=0; j<M; j++)
                    Pmat (i ,j) = _tt_cov.at (i+min_i).at (j+min_i);
#if VERBOSE == 2
            std::cout << "\n P = \n" << Pmat;
#endif

            // construct F matrix
            Eigen::MatrixXd F (M, M);
            F.setZero ();
            for (int i=0; i<M; i++)
            {
                // if no time has passed, we shouldn't *actually* be getting here!
                if (tt_delta == 0) F (i, i) = 1;
                // if current estimate is LESS THAN ZERO, prediction should also be 0
                else if (Bmat(i) <= 1e-6) F (i, i) = 0;
                else F (i, i) = fmin(1, fmax(0, 1 - (double)tt_delta / Bmat (i)));
            }
#if VERBOSE == 2
            std::cout << "\n F = \n" << F.format (decimalMatFmt);
#endif
            
            // Identity matrix
            Eigen::MatrixXd I (M, M);
            I = Eigen::MatrixXd::Identity (M, M);

            // System noise matrix (variability/second)
            Eigen::MatrixXd Q = I;
            for (int i=0; i<M; i++) Q (i, i) = tt_delta;
#if VERBOSE == 2
            std::cout << "\n Q = \n" << Q.format (decimalMatFmt);
#endif
            

            // Predict Bhat
            Eigen::VectorXd Bhat (M);
            Bhat = F * Bmat;
#if VERBOSE == 2
            std::cout << "\n Bhat = \n" << Bhat;
#endif

            // Predict Phat 
            Eigen::MatrixXd Phat (M, M);
            Phat = F * Pmat * F.transpose () + Q;
#if VERBOSE == 2
            std::cout << "\n Phat = \n" << Phat;
#endif

            // --- Update equations
            // Convert Z to (column) vector
            Eigen::VectorXd Z (M);
            for (int i=0; i<M; i++) Z (i) = tt_obs.at (i+min_i);
            // std::cout << "\n Z = \n" << Z;

            // Convert R to matrix
            Eigen::MatrixXd R (M, M);
            for (int i=0; i<M; i++)
                for (int j=0; j<M; j++)
                    R (i, j) = tt_r.at (i+min_i).at (j+min_i);
            // std::cout << "\n R = \n" << R;

            
            // Measurement matrix
            Eigen::MatrixXd H = I;
            // lower diagonal is -1's
            for (int i=1; i<M; i++)
                    H (i, i-1) = -1;
#if VERBOSE == 2
            std::cout << "\n Measurement matrix H = \n" << H;
#endif

            // Innovation residual
            Eigen::VectorXd resid_y = Z - H * Bhat;
#if VERBOSE == 2
            std::cout << "\n HBhat = \n" << (H * Bhat);
            std::cout << "\n yresid = \n" << resid_y;
#endif

            // Innovation covariance
            Eigen::MatrixXd S = R + H * Phat * H.transpose ();
            // std::cout << "\n S = \n" << S;

            // Kalman gain
            Eigen::MatrixXd K = Phat * H.transpose () * S.inverse ();
#if VERBOSE == 2
            std::cout << "\n K = \n" << K.format (decimalMatFmt);
#endif

            // Updated state estimate
            Bhat = Bhat + K * resid_y;
#if VERBOSE == 2
            std::cout << "\n Bhat = \n" << Bhat;
#endif

            // Update state covariance
            Phat = (I - K * H) * Phat * (I - K * H).transpose () +
                K * R * K.transpose ();
#if VERBOSE == 2
            std::cout << "\n Phat = \n" << Phat;
#endif

            // Save state
            int Mx = stops.size ();
            for (int i=0; i<Mx; i++)
            {
                _tt_state.at (i) = (i < min_i) ? 0 : Bhat (i - min_i);
            }
#if VERBOSE == 2
            std::cout << "\n\n New state vector: [";
            for (auto z : _tt_state) std::cout << " " << z << " ";
            std::cout << "]";
#endif

            for (int i=0; i<Mx; i++)
            {
                for (int j=0; j<Mx; j++)
                {
                    _tt_cov.at (i).at (j) = (i < min_i || j < min_i) ? 0 : Phat (i-min_i, j-min_i);
                }
            }
#if VERBOSE == 2
            std::cout << "]\n with covariance matrix: ";
            for (auto br : _tt_cov)
            {
                std::cout << "\n  ";
                for (auto bc : br) std::cout << std::setw(9) << std::round (bc) << "  ";
            }
            std::cout << "\n\n --- the correlation matrix for this looks like ... ";
            for (int i=0; i<Mx; i++)
            {
                std::cout << "\n  ";
                for (int j=0; j<Mx; j++)
                {
                    if (j == i) 
                    {
                        std::cout << std::setw (9) << "1  ";
                    }
                    else if (_tt_cov.at (i).at (j) == 0)
                    {
                        std::cout << std::setw (9) << "0  ";
                    }
                    else
                    {
                        std::cout << std::setw(9)
                            << (_tt_cov.at (i).at (j) / 
                                pow (_tt_cov.at (i).at (i), 0.5) / 
                                pow (_tt_cov.at (j).at (j), 0.5)) << "  ";
                    }
                }
            }
            std::cout << "\n\n";
#endif
        

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

        int curstop = _current_stop;

        double tsum, tvar;
        // calculate the timestamp of the scheduled start time 
        uint64_t t0 = _timestamp - (Time (_timestamp) - this->trip_start_time ());
        for (int i=0; i<M; ++i)
        {
            etas.at (i).stop_id = stops.at (i).stop->stop_id ();
            etas.at (i).estimate = 0;
            
            // it should be of the time the vehicle ETAs were last updated 
            // (which should always be the same time as the vehicle's timestamp ...)
            // tsum = 0.0;
            // tvar = 0.0;
            // for (int j=0; j<i; j++)
            // {
            //     tsum += _tt_state.at (j);
            //     for (int k=0; k<i; k++)
            //     {
            //         tvar += _tt_cov.at (j).at (k);
            //     }
            // }


            double tsum, tvar;
            tsum += _tt_state.at (i);
            if (tsum <= 0) continue;

            tvar = 0;
            for (int j=0;j<=i;j++) for (int k=0; k<=i; k++) tvar += _tt_cov.at (j).at (k);
            
            if (i <= curstop) continue;
            etas.at (i).estimate = t0 + tsum;
            // for now just the 95% credible interval
            etas.at (i).quantiles.emplace_back (0.025, t0 + (tsum - 2 * pow(tvar, 0.5)));
            etas.at (i).quantiles.emplace_back (0.975, t0 + (tsum + 2 * pow(tvar, 0.5)));
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

    }

}
