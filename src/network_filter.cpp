#include "network_filter.h"

namespace Gtfs {

    std::pair<double,double> Segment::predict (int delta)
    {
        // use current estimate and (historical) prior to predict future state
        double xhat, Phat;
        xhat = _travel_time;
        Phat = _uncertainty + delta * 2.0 / 60;

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

} // end Gtfs 
