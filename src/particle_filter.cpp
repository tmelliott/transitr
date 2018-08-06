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
            p->travel (_delta);
        }
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
        std::cout << "\nsumlh = " << sumlh << " from "
            << _state.size () << " particles";

        // compute particle (cumulative) weights ll - log(sum(likelihood)))
        

        // resample u \in [0, 1)
        
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
        return 100 * d / dmax + 0.5;
    }



    Particle::Particle (double d, double s, Vehicle* v)
    {
        vehicle = v;
        distance = d;
        speed = s;
    }

    double Particle::get_distance ()
    {
        return distance;
    }

    double Particle::get_speed ()
    {
        return speed;
    }

    void Particle::travel (unsigned delta)
    {
        // do the particle physics
        double Dmax = vehicle->trip ()->shape ()->path ().back ().distance;
        while (distance < Dmax && delta > 0)
        {
            distance += speed;
            delta--;
        }
    }

    double 
    Particle::calculate_likelihood (latlng& y, 
                                    std::vector<ShapePt>* path, 
                                    double sigma)
    {
        log_likelihood = - log (sigma);
        return exp (log_likelihood);
    }

}; // namespace Gtfs
