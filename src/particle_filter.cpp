/**
 * Define the particle filter functions.
 */
#include "particle_filter.h"
#include "timing.h"

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
#if WRITE_PARTICLES
            std::ostringstream fname;
            fname << "history/vehicle_" << _vehicle_id << ".csv";
            std::ofstream fout;
            fout.open (fname.str ().c_str (), std::ofstream::app);
            double d (p->get_distance ());
            latlng ppos (_trip->shape ()->coordinates_of (d));
            fout << _timestamp << ","
                << _trip->trip_id () << ","
                << _position.latitude << "," << _position.longitude << ","
                << p->get_distance () << ","
                << p->get_speed () << ","
                << std::setprecision(15)
                << ppos.latitude << "," << ppos.longitude << "\n";
            fout.close ();

            std::ostringstream fname2;
            fname2 << "history/vehicle_" << _vehicle_id << "_particles.csv";
            fout.open (fname2.str ().c_str (), std::ofstream::app);
            fout << _timestamp;
            for (auto ati = p->get_arrival_times ().begin (); ati != p->get_arrival_times ().end (); ++ati)
            {
                fout << "," << *ati;
            }
            fout << "\n";
            fout.close ();
#endif
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

    void Vehicle::predict_etas (RNG& rng)
    {
        if (!valid ()) return;

#if VERBOSE == 2
        Timer timer;
        std::cout << "- vehicle " << _vehicle_id << " - predicting etas";
#endif
        for (auto p = _state.begin (); p != _state.end (); ++p) p->predict_etas (rng);
#if VERBOSE == 2
        std::cout << " (" << timer.cpu_seconds () << "ms)\n";
#endif
    }

    etavector Vehicle::get_etas ()
    {
        auto stops = _trip->stops ();
        int M (stops.size ());
        etavector etas;
        etas.resize (M);
        for (int i=0; i<M; ++i)
        {
            // need to center each particle's arrival time
            int tarr = 0;
            int ni = 0;
            for (auto p = _state.begin (); p != _state.end (); ++p)
            {
                if (p->get_arrival_time (i) > 0)
                {
                    tarr += (p->get_arrival_time (i) - _timestamp);
                    ni++;
                }
            }
            etas.at (i).stop_id = stops.at (i).stop->stop_id ();
            if (ni == 0) continue;
            tarr /= ni;
            etas.at (i).estimate = _timestamp + tarr;
        }
        return etas;
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
        if (distance >= Dmax) 
        {
            distance = Dmax;
            complete = true;
            return;
        }
        
        // get STOPS
        std::vector<StopTime>* stops;
        if (!vehicle || !vehicle->trip ()) return;
        stops = &(vehicle->trip ()->stops ());
        int M (stops->size ());
        unsigned int m (find_stop_index (distance, stops));
        if (m == M-1) 
        {
            distance = Dmax;
            complete = true;
            return;
        }

        double next_stop_d = stops->at (m).distance;
        
        // get SEGMENTS
        

        while (distance < Dmax && delta > 0)
        {
            // add noise to speed
            double speed_prop = -1;
            while (speed_prop < 0 || speed_prop > 30)
            {
                speed_prop = speed + rng.rnorm () * (5 / 30);
            }
            speed = speed_prop;
            if (distance + speed >= next_stop_d)
            {
                // about to reach a stop ... slow? stop? just drive past?
                at.at (m) = vehicle->timestamp () - delta - 1;
                if (rng.runif () < 0.5)
                {
                    // stop dwell time ~ Exp(tau = 10)
                    double gamma = 6;
                    double tau = 10;
                    double dwell = gamma - tau * log (rng.runif ());
                    delta = fmax(0, delta - dwell);
                    distance = next_stop_d;
                    if (m == M-1) break;
                    m++;
                    next_stop_d = stops->at (m).distance;
                    continue;
                }
            }
            distance += speed;
            delta--;
        }
    }

    void Particle::predict_etas (RNG& rng)
    {
        if (complete) return;

        // get STOPS
        std::vector<StopTime>* stops;
        if (!vehicle || !vehicle->trip ()) return;
        stops = &(vehicle->trip ()->stops ());
        int M (stops->size ());
        unsigned int m (find_stop_index (distance, stops));
        if (m == M-1) 
        {
            return;
        }

        double dcur = distance;
        uint64_t t0 = vehicle->timestamp ();
        while (m < M)
        {
            double dnext = stops->at (m).distance;
            at.at (m) = t0 + (dnext - dcur) / speed;
            dcur = dnext;
            // and add some dwell time
            t0 = at.at (m);
            m++;
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
