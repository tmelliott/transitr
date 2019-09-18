#include "network_filter.h"

namespace Gtfs {

    void Segment::update (par* params, Gtfs* gtfs)
    {
        if (!loaded) load ();
        // set up the initial state of the segment
        _travel_time (0) = _length / 10.0;
        _travel_time (1) = 0.0;
        _uncertainty = Eigen::Matrix2d::Zero ();
        min_tt = _length / max_speed;

        // also specify the "between vehicle" variabilty
        // _state_var = pow (log (10 * min_tt), 2);

        if (_system_noise == 0)
            _system_noise = params->nw_system_noise;

        // best estimate
        // log(phi_ell) = theta0 + theta1 * log(min_tt)
        if (_state_var == 0)
            _state_var = exp (-1.2 + 1.2 * log (min_tt));
        
        _measurement_error = params->nw_measurement_error;

#if VERBOSE > 1
        // std::cout << "\n - using parameters [q, phi, e] = ["
        //     << _system_noise
        //     << ", " << _state_var 
        //     << ", " << _measurement_error << "]\n\n";
#endif
    }

    std::pair<Eigen::Vector2d, Eigen::Matrix2d> Segment::predict (int delta)
    {
        if (!loaded) load ();
        // use current estimate and (historical) prior to predict future state
        auto xhat = _travel_time;
        auto Phat = _uncertainty;

        if (_uncertainty (0, 0) > 0)
        {
            Phat (0, 0) += pow (_system_noise * delta, 2);
            // Eigen::Matrix2d F, Q;
            // switch (_model_type) 
            // {
            //     case 0:
            //         F << 1, 0, 0, 0;
            //         Q << pow (_system_noise * delta, 2), 0, 0, 0;
            //         break;
            //     case 1:
            //         // I think this is wronggg? squares etc.
            //         F << 1, delta, 0, 1;
            //         Q << 0, 0, 0, _system_noise;
            //         break;
            // }
            // Phat = F * Phat * F.transpose () + Q;
        }
        else
        {
            Phat = Eigen::Matrix2d::Zero ();
            switch (_model_type) 
            {
                case 0:
                    Phat (0, 0) = 1000.0;
                    break;
                case 1:
                    Phat (0, 0) = 100.0; // fix this
                    Phat (1, 1) = 100.0;
                    break;
            }
        }

        return std::make_pair (xhat, Phat);
    }

    std::pair<Eigen::Vector2d, Eigen::Matrix2d> Segment::predict (uint64_t t)
    {
        if (t <= _timestamp) return predict (0);
        return predict ((int)(t - _timestamp));
    }

    void Segment::update (uint64_t now)
    {
        if (!loaded) load ();
        if (_data.size () == 0) return;

        // first, predict the future state ...
        if (_timestamp == 0)
        {
            _timestamp = now;
        }
        
#if VERBOSE > 0
        std::cout << "\n\n + Segment " << _segment_id 
            << "\n  => X = "<< _travel_time.format (tColVec) 
            << " and P = " << _uncertainty.format(inlineMat);
#endif

        Eigen::Vector2d xhat = _travel_time;
        Eigen::Matrix2d Phat = _uncertainty;
        if (Phat (0, 0) == 0) Phat (0, 0) = 1000;
        // auto res = predict (now - _timestamp);
        // xhat = res.first;
        // Phat = res.second;
        int delta = now - _timestamp;
        if (delta > 0) 
        {
            Phat (0, 0) += pow(delta * _system_noise, 2);
        }

        double B, b, U, u, Z, z;
        B = Phat (0, 0);
        b = xhat (0);
        U = 1 / B;
        u = b / B;
        Z = 0;
        z = 0;
        double err;
        std::pair<int, double>* dj;
        for (int j=0; j<_data.size (); j++)
        {
            dj = &(_data.at (j));
            err = fmax (_measurement_error, dj->second) + pow (_state_var, 2);
            Z += 1 / err;
            z += dj->first / err;
        }
        U += Z;
        u += z;

        B = 1 / U;
        b = u / U;

        Phat (0, 0) = B;
        xhat (0) = b;

// #if VERBOSE > 0
//         std::cout
//             << "\n  => X = "<< xhat.format (tColVec) 
//             << " and P = " << Phat.format(inlineMat);
// #endif

        
//         // then update with observations
//         {
//             // transform to information 
//             Eigen::Matrix2d Z = Eigen::Matrix2d::Zero ();
//             Eigen::Vector2d z = Eigen::Vector2d::Zero ();
//             switch (_model_type) 
//             {
//                 case 0:
//                     // only worrying about single dimension
//                     Z (0, 0) = pow (Phat (0, 0), -1);
//                     z (0, 0) = xhat (0, 0) / Phat (0, 0);
//                     break;
//                 case 1:
//                     Z = Phat.inverse ();
//                     z = Z * xhat;
//                     break;
//             }

// #if VERBOSE > 0
//             std::cout << "\n  => U ="
//                 << Z.format (inlineMat) << " and i = " << z.format (tColVec);
// #endif

//             // this is the same regardless of the model being used
//             Eigen::MatrixXd H (1, 2);
//             H << 1, 0;

//             double err = _measurement_error;
//             Eigen::Matrix2d I = Eigen::Matrix2d::Zero ();
//             Eigen::Vector2d i = Eigen::Vector2d::Zero ();
//             std::pair<int, double>* dj;
//             double errj;
//             for (int j=0; j<_data.size (); j++)
//             {
//                 dj = &(_data.at (j));
//                 // "err" is the sort-of minimum error we expect 
//                 // for segment's travel time observations
//                 errj = fmax (err, pow (dj->second, 2)) + pow (_state_var, 2);
//                 // errj = fmin (err, fmax (dj->first, dj->second)) + _state_var;
// #if VERBOSE > 0
//                 std::cout 
//                     << "\n  => Z = " << dj->first 
//                     << ", R = (" << dj->second << " \u2227 " << err
//                     << ") + " << _state_var
//                     << " -> err = " << errj;
// #endif
//                 I += H.transpose () * (1 / errj) * H;
//                 i += H.transpose () * (dj->first / errj);
//             }

// #if VERBOSE > 0
//             std::cout << "\n  => I ="
//                 << I.format (inlineMat) << " and i = " << i.format (tColVec);
// #endif

//             Z += I;
//             z += i;

// #if VERBOSE > 0
//             std::cout << "\n  => U ="
//                 << Z.format (inlineMat) << " and i = " << z.format (tColVec);
// #endif

//             // reverse transform information to travel time state space
//             switch (_model_type) 
//             {
//                 case 0:
//                     Phat (0, 0) = pow (Z (0, 0), -1);
//                     xhat (0, 0) = z (0, 0) / Z (0, 0);
//                     break;
//                 case 1:
//                     Phat = Z.inverse ();
//                     xhat = Phat * z;
//                     break;
//             }

// #if VERBOSE > 0
//             std::cout << "\n  => X = "
//                 << xhat.format (tColVec) 
//                 << " and P = " << Phat.format(inlineMat);
// #endif
//         }

        _travel_time = xhat;
        _uncertainty = Phat;
        _timestamp = now;
        _data.clear ();
    }

    double Segment::get_speed ()
    {
        if (!loaded) load ();
        if (_travel_time (0) > 0)
        {
            return _length / _travel_time (0);
        }
        return 0.0;
    }
    double Segment::get_speed (int delta)
    {
        if (!loaded) load ();
        if (_travel_time (0) > 0)
        {
            auto x = predict (delta);
            return _length / x.first (0);
        }
        return 0.0;
    }
    int Segment::sample_travel_time (RNG& rng)
    {
        sample_travel_time (rng, 0);
    }
    int Segment::sample_travel_time (RNG& rng, int delta)
    {
        if (!loaded) load ();
        if ( _travel_time (0) == 0) 
        {
            return _length / (rng.rnorm () * 25.0 + 5.0); // [5, 30m/s]
        }

        auto x = predict (delta);

        // truncated normal distribution
        double tt (-1.0);
        int tries (100);
        while (tt < min_tt && tries > 0) {
            tt = rng.rnorm () * (x.second (0, 0) + _state_var) + x.first (0);
            tries--;
        }
        return round (fmax (0.0, tt));
    }
    double Segment::sample_speed (RNG& rng)
    {
        sample_speed (rng, 0);
    }
    double Segment::sample_speed (RNG& rng, int delta)
    {
        int x = sample_travel_time (rng, delta);
        x = std::max(5, x);
        return _length / x;
    }

} // end Gtfs 
