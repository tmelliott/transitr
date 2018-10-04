#include "network_filter.h"

namespace Gtfs {

    std::pair<double,double> Segment::predict (int delta)
    {
        // use current estimate and (historical) prior to predict future state
        double xhat, Phat;
        xhat = _travel_time;
        Phat = _uncertainty + (double) delta * 2.0 / 60.0 / 30.0;

        return std::make_pair (xhat, Phat);
    }

    void Segment::update (uint64_t now)
    {
        if (_data.size () == 0) return;

        // first, predict the future state ...
        double xhat, Phat;
        if (_uncertainty == 0)
        {
            xhat = _length / 10.0;
            Phat = 1000;
        }
        else
        {
            auto res = predict (now - _timestamp);
            xhat = res.first;
            Phat = res.second;
        }

        // std::cout << "\n + Segment " << _segment_id << " >> "
        //     << "[" << xhat << ", " << Phat << "] >> ";

        
        // then update with observations
        {
            // transform to information 
            double Z = pow(Phat, -1);
            double z = xhat * Z;

            double I = std::accumulate (_data.begin (), _data.end (), 0.0, 
                                        [](double a, std::pair<int, double>& b) {
                                            return a + pow (b.second, -1);
                                        });
            // std::cout << "(";
            double i = std::accumulate (_data.begin (), _data.end (), 0.0, 
                                        [](double a, std::pair<int, double>& b) {
                                            // std::cout << b.first << ", " << b.second << " ; ";
                                            return a + (double) b.first * pow (b.second, -1);
                                        });
            // std::cout << ") >> ";

            Z += I;
            z += i;

            // reverse transform information to travel time state space
            Phat = pow (Z, -1);
            xhat = z * Phat;
        }

        // std::cout << "[" << xhat << ", " << Phat << "]";

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
    double Segment::get_speed (int delta)
    {
        if (_uncertainty > 0 && _travel_time > 0)
        {
            auto x = predict (delta);
            return _length / x.first;
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
