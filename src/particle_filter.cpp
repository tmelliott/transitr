/**
 * Define the particle filter functions.
 */
#include "particle_filter.h"

namespace Gtfs {

    void Vehicle::initialize ()
    {
        _state.clear ();
        _state.reserve (_N);
        double dmax = _trip->shape ()->path ().back ().distance;
        // std::cout << "\n + Initialize " << _vehicle_id << ": ";
        for (int i=0; i<_N; ++i)
        {
            // initialize each particle with a state X(i)
            _state.emplace_back (10.0, 1.0, this);
            // std::cout << " [" << _state.back ().get_distance () << "]";
        }

        _newtrip = false;
    }

    void Vehicle::mutate ()
    {
        if (_newtrip)
        {
            initialize ();
            return;
        }

        // There probably need to be a bunch of checks here ...
        
        // do the transition ("mutation")
        for (auto p = _state.begin (); p != _state.end (); ++p)
        {
            p->travel (_delta);
        }
    }

    void Vehicle::select ()
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
