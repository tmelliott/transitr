#include "network_filter.h"

namespace Gtfs {

    std::pair<double,double> Segment::predict (int delta)
    {
        // use current estimate and (historical) prior to predict future state
        double xhat, Phat;
        xhat = _travel_time;
        // Phat = _uncertainty + delta * 2.0 / 60;
        Phat = _uncertainty + 2.0;

        return std::make_pair (xhat, Phat);
    }

    void Segment::update (uint64_t now)
    {
        if (_data.size () == 0) return;

        // first, predict the future state ...
        double xhat, Phat;
        if (_uncertainty == 0)
        {
            xhat = 15;
            Phat = 10;
        }
        else
        {
            auto res = predict (now - _timestamp);
            xhat = res.first;
            Phat = res.second;
        }
        
        // then update with observations
        {
            double y = std::accumulate (_data.begin (), _data.end (), 0.0);
            y /= (double) _data.size ();
            double R = 5.0;

            double z = y - xhat;
            double S = R + Phat;
            double K = Phat * pow(S, -1);

            xhat += K * z;
            Phat *= 1 - K;
        }

        _travel_time = xhat;
        _uncertainty = Phat;
        _timestamp = now;
        _data.clear ();
    }

    double Segment::get_speed ()
    {
        if (_uncertainty > 0 && _travel_time > 0)
        {
            return _length / _travel_time;
        }
        return 0.0;
    }
    int Segment::sample_travel_time (RNG& rng)
    {
        if (_uncertainty > 0 && _travel_time > 0)
        {
            double x = rng.rnorm () * _uncertainty + _travel_time;
            return round (fmax (0.0, x));
        }
        return 0.0;
    }
    double Segment::sample_speed (RNG& rng)
    {
        int x = sample_travel_time (rng);
        if (x == 0) return 0.0;
        return _length / x;
    }

} // end Gtfs 
