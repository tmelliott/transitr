/**
 * Define the particle filter functions.
 */
#include "particle_filter.h"

namespace Gtfs {

    void Vehicle::initialize (RNG& rng)
    {
        _state.clear ();
        _state.reserve (_N);
        double dmax = _trip->shape ()->path ().back ().distance;
        // std::cout << "\n + Initialize " << _vehicle_id << " [N=" << _N << "]: ";
        for (int i=0; i<_N; ++i)
        {
            // initialize each particle with a state X(i) - distance, speed
            _state.emplace_back (rng.runif () * dmax,
                                 rng.runif () * 30, 
                                 this);
            // std::cout << " [" << _state.back ().get_distance () << ", " 
            //     << _state.back ().get_speed () << "]";
        }

        _newtrip = false;
    }

    void Vehicle::mutate (RNG& rng)
    {
        if (_newtrip)
        {
            initialize (rng);
            return;
        }

        if (!valid () || _delta == 0) return;

        // There probably need to be a bunch of checks here ...
        
        // do the transition ("mutation")
        for (auto p = _state.begin (); p != _state.end (); ++p)
        {
            p->travel (_delta, rng);
        }

        // update
        select (rng);
    }

    void Vehicle::select (RNG& rng)
    {
        // calculate loglikelihood of particles
        double sumlh = 0.0;
        std::vector<ShapePt>* path = &(_trip->shape ()->path ());
        for (auto p = _state.begin (); p != _state.end (); ++p)
        {
            sumlh += p->calculate_likelihood (_position, path, _gpserror);
        }

        // compute particle (cumulative) weights ll - log(sum(likelihood)))
        std::vector<double> wt;
        wt.reserve (_N+1);
        wt.push_back (0);
        double lsumlh (log (sumlh));
        for (auto p = _state.begin (); p != _state.end (); ++p)
        {
            wt.push_back (wt.back () + exp (p->get_ll () - lsumlh));
        }

        // some condition for resampling (1/wt^2 < N? or something ...)


        if (wt.back () < 0.99)
        {
            // bad bad bad
            _state.clear ();
            _newtrip = true;
            return;
        }

        // resample u \in [0, 1)
        double u;
        unsigned int j;
        std::vector<Particle> _newstate;
        _newstate.reserve (_N);
        for (int i=0; i<_N; ++i)
        {
            u = rng.runif ();
            j = 1;
            while (wt[j] <= u && j < _N) j++;
            _newstate.emplace_back (_state[j-1]);
        }

        // std::cout << "\n + vehicle " << _vehicle_id << ":\n";
        // for (auto p = _state.begin (); p != _state.end (); ++p)
        //     std::cout << "[" << p->get_distance () << ", " << p->get_speed () << "] ; ";
        _state.clear ();
        _state = std::move(_newstate);
        // std::cout << "\n";
        // for (auto p = _state.begin (); p != _state.end (); ++p)
        //     std::cout << "[" << p->get_distance () << ", " << p->get_speed () << "] ; ";
    }

    double Vehicle::distance ()
    {
        double distance = 0.0;
        for (auto p = _state.begin (); p != _state.end (); ++p) distance += p->get_distance ();
        return distance / (double)_state.size ();
    }

    double Vehicle::speed ()
    {
        double speed = 0.0;
        for (auto p = _state.begin (); p != _state.end (); ++p) speed += p->get_speed ();
        return speed / (double)_state.size ();
    }

    int Vehicle::progress ()
    {
        double d = distance ();
        double dmax = _trip->shape ()->path ().back ().distance;
        return (100 * d / dmax) + 0.5;
    }



    void Particle::travel (unsigned delta, RNG& rng)
    {
        // do the particle physics
        double Dmax = vehicle->trip ()->shape ()->path ().back ().distance;
        while (distance < Dmax && delta > 0)
        {
            double speed_prop = -1;
            while (speed_prop < 0 || speed_prop > 30)
            {
                speed_prop = speed + rng.rnorm () * 2.0;
            }
            speed = speed_prop;
            distance += speed;
            delta--;
        }
    }

    double 
    Particle::calculate_likelihood (latlng& y, 
                                    std::vector<ShapePt>* path, 
                                    double sigma)
    {
        latlng ppos = vehicle->trip ()->shape ()->coordinates_of (distance);
        // (log) distance between points
        double ld = log (distanceEarth (ppos, vehicle->position ()));
        // (log) (d/sigma)^2 ~ Chi2(2) ~ Exp(2)
        double lX2 = 2 * (ld - log (sigma));
        // log pdf of lX2 ~ Exp(2)
        log_likelihood = log (0.5) - 0.5 * exp(lX2);
        return exp (log_likelihood);
    }

}; // namespace Gtfs
