#include "network_filter.h"

namespace Gtfs {

    void Segment::update (par* params, Gtfs* gtfs)
    {
        if (!loaded) load ();
        // set up the initial state of the segment
        _speed = _prior_speed;
        _uncertainty = _prior_speed_var;

        if (_system_noise == 0)
            _system_noise = params->nw_system_noise;

        // also specify the "between vehicle" variabilty
        if (_state_var == 0)
            _state_var = 5. / pow (3.6, 2);

        _measurement_error = params->nw_measurement_error;
    }

    std::pair<double, double> Segment::predict (int delta)
    {
        if (!loaded) load ();
        // use current estimate and (historical) prior to predict future state
        double xhat = _speed;
        double Phat = _uncertainty;

        if (_uncertainty > 0)
        {
            // Check this ... (chapters 4/5):
            Phat += pow (_system_noise * delta, 2);
            if (delta > 30 && _prior_speed_var > 0 &&
                Phat > 2 * _prior_speed_var)
                Phat = 2 * _prior_speed_var;
        }
        else
        {
            Phat = 0.;
            switch (_model_type)
            {
                case 0:
                    Phat = 50. / 4.; // approx 50~km/h/s
                    break;
                case 1:
                    Phat = 50. / 4.; // fix this
                    break;
            }
        }

        return std::make_pair (xhat, Phat);
    }

    std::pair<double, double> Segment::predict (uint64_t t)
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

        bool log (false);
#if VERBOSE > 0
        log = true;
        std::cout << "\n\n + Segment " << std::setw (6) << _segment_id
            << ", delta = " << std::setw (4) << (now - _timestamp)
            << " -> X = " << std::setw (5) << round (_speed)
            << ", P = " << std::setw (6) << round (_uncertainty);
#endif

        double xhat = _speed;
        double Phat = _uncertainty;
        if (Phat == 0.) Phat = 5.;
        int delta = now - _timestamp;
        if (delta > 0.)
        {
            Phat += pow (delta * _system_noise, 2.);
        }

        if (log)
            std::cout
                << " => X = " << std::setw (5) << round (_speed)
                << ", P = " << std::setw (6) << round (_uncertainty);

        double B, b, U, u, Z, z;
        B = Phat;
        b = xhat;
        U = 1 / B;
        u = b / B;
        Z = 0;
        z = 0;
        double err;
        std::pair<double, double>* dj;
        if (log) std::cout << "; Z: [";
        for (int j=0; j<_data.size (); j++)
        {
            dj = &(_data.at (j));
            err = fmax (_measurement_error, dj->second) + pow (_state_var, 2);
            if (log)
                std::cout << "{"
                    << round (dj->first) << ", " << round (err)
                    << "}";
            Z += 1 / err;
            z += dj->first / err;
        }
        if (log) std::cout << "]";
        U += Z;
        u += z;

        B = 1 / U;
        b = u / U;

        if (log)
            std::cout
                << " => X = " << std::setw (5) << round (b)
                << ", P = " << std::setw (6) << round (B);

        Phat = B;
        xhat = b;

        _speed = xhat;
        _uncertainty = Phat;
        _timestamp = now;
        _data.clear ();
    }

    double Segment::speed (int delta)
    {
        if (!loaded) load ();
        if (_speed > 0.)
        {
            auto x = predict (delta);
            return x.first;
        }
        return 0.0;
    }

    double Segment::sample_speed (RNG& rng)
    {
        sample_speed (rng, 0);
    }
    double Segment::sample_speed (RNG& rng, int delta)
    {
        int nmax = 100;
        double s = 0.0;
        auto x = predict (delta); // [speed, uncertainty]
        while (s < 0.5 || s > _max_speed)
        {
            s = rng.rnorm () * x.second + x.first;
            if (nmax-- == 0)
            {
                return rng.runif () * (_max_speed - 10.) + 5.;
            }
        }
        return s;
    }

    int Segment::sample_travel_time (RNG& rng)
    {
        sample_travel_time (rng, 0.);
    }
    int Segment::sample_travel_time (RNG& rng, int delta)
    {
        if (!loaded) load ();
        double speed = sample_speed (rng, delta);
        return round (_length / speed);
    }

} // end Gtfs
