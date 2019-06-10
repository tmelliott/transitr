#include "network_filter.h"

namespace Gtfs {

    void Segment::update (par* params, Gtfs* gtfs)
    {
        if (!loaded) load ();
        // set up the initial state of the segment
        _system_noise = params->nw_system_noise;
        _measurement_error = params->nw_measurement_error;

        // fetch prior from the database
        // ...
        
        _travel_time (0) = _length / 10.0;
        _travel_time (1) = 0.0;
        _uncertainty = Eigen::Matrix2d::Zero ();
        min_tt = _length / max_speed;
    }

    std::pair<Eigen::Vector2d, Eigen::Matrix2d> Segment::predict (int delta)
    {
        if (!loaded) load ();
        // use current estimate and (historical) prior to predict future state
        auto xhat = _travel_time;
        auto Phat = _uncertainty;

        if (_uncertainty (0, 0) > 0)
        {
            Eigen::Matrix2d F, Q;
            F << 1, delta, 0, 1;
            Q << 0, 0, 0, _system_noise;
            Phat = F * Phat * F.transpose () + Q;
        }
        else
        {
            Phat = Eigen::Matrix2d::Zero ();
            Phat (0, 0) = _travel_time (0);
            Phat (1, 1) = 100.0;
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

        Eigen::Vector2d xhat;
        Eigen::Matrix2d Phat;
        auto res = predict (now - _timestamp);
        xhat = res.first;
        Phat = res.second;


#if VERBOSE > 0
        std::cout
            << "\n  => X = "<< xhat.format (tColVec) 
            << " and P = " << Phat.format(inlineMat);
#endif

        
        // then update with observations
        {
            // transform to information 
            Eigen::Matrix2d Z = Phat.inverse ();
            Eigen::Vector2d z = Z * xhat;

            Eigen::MatrixXd H (1, 2);
            H << 1, 0;

            double err = _measurement_error;
            Eigen::Matrix2d I = Eigen::Matrix2d::Zero ();
            Eigen::Vector2d i = Eigen::Vector2d::Zero ();
            std::pair<int, double>* dj;
            double errj;
            for (int j=0; j<_data.size (); j++)
            {
                dj = &(_data.at (j));
                errj = fmin (err, fmax (dj->first, dj->second));
#if VERBOSE > 0
                std::cout << "\n  => Z = "
                    << dj->first << ", R = " << dj->second
                    << " -> err = " << errj;
#endif
                I += H.transpose () * (1 / errj) * H;
                i += H.transpose () * (dj->first / errj);
            }

#if VERBOSE > 0
            std::cout << "\n  => I ="
                << I.format (inlineMat) << " and i = " << i.format (tColVec);
#endif

            Z += I;
            z += i;

            // reverse transform information to travel time state space
            Phat = Z.inverse ();
            xhat = Phat * z;
#if VERBOSE > 0
            std::cout << "\n  => X = "
                << xhat.format (tColVec) 
                << " and P = " << Phat.format(inlineMat);
#endif
        }

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
            tt = rng.rnorm () * x.second (0, 0) + x.first (0);
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
