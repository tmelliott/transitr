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
        
        _travel_time = _length / 10.0;
        _uncertainty = 0.0;
        // _uncertainty = _travel_time; // m/s
        min_tt = _length / max_speed;
    }

    std::pair<double,double> Segment::predict (int delta)
    {
        if (!loaded) load ();
        // use current estimate and (historical) prior to predict future state
        double xhat, Phat;
        xhat = _travel_time;

        if (_uncertainty > 0)
        {
            Phat = _uncertainty + delta * _system_noise;
        }
        else
        {
            Phat = _travel_time;
        }

        return std::make_pair (xhat, Phat);
    }

    std::pair<double,double> Segment::predict (uint64_t t)
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
        
        double xhat, Phat;
        auto res = predict (now - _timestamp);
        xhat = res.first;
        Phat = res.second;

#if VERBOSE > 0
        std::cout << "\n + Segment " << _segment_id << " >> "
            << "[" << xhat << ", " << Phat << "] >> ";
#endif

        
        // then update with observations
        {
            // transform to information 
            double Z = pow(Phat, -1);
            double z = xhat * Z;

#if VERBOSE > 0
            std::cout << "(";
#endif
            double err = _measurement_error;
            double I = std::accumulate (
                _data.begin (), 
                _data.end (), 
                0.0, 
                [&err](double a, std::pair<int, double>& b) 
                {
                    return a + pow (fmin (err, fmax (b.second, b.first)), -1);
                }
            );
            double i = std::accumulate (
                _data.begin (), 
                _data.end (), 
                0.0, 
                [&err](double a, std::pair<int, double>& b) 
                {
#if VERBOSE > 0
                    std::cout << b.first << ", " << b.second << " ; ";
#endif
                    return a + (double) b.first * 
                        pow (fmin (err, fmax (b.second, b.first)), -1);
                }
            );
#if VERBOSE > 0
            std::cout << ") >> { "
                << I << ", " << i << " } >> ";
#endif

            Z += I;
            z += i;

            // reverse transform information to travel time state space
            Phat = pow (Z, -1);
            xhat = z * Phat;
        }

#if VERBOSE > 0
        std::cout << "[" << xhat << ", " << Phat << "]";
#endif

        _travel_time = xhat;
        _uncertainty = Phat;
        _timestamp = now;
        _data.clear ();
    }

    double Segment::get_speed ()
    {
        if (!loaded) load ();
        if (_travel_time > 0)
        {
            return _length / _travel_time;
        }
        return 0.0;
    }
    double Segment::get_speed (int delta)
    {
        if (!loaded) load ();
        if (_travel_time > 0)
        {
            auto x = predict (delta);
            return _length / x.first;
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
        if ( _travel_time == 0) 
        {
            return _length / (rng.rnorm () * 25.0 + 5.0); // [5, 30m/s]
        }

        auto x = predict (delta);

        // truncated normal distribution
        double tt (-1.0);
        int tries (100);
        while (tt < min_tt && tries > 0) {
            tt = rng.rnorm () * x.second + x.first;
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
