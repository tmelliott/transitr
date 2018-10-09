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
            _state.emplace_back ((double)i / (double)(_N - 1) * dmax, //rng.runif () * dmax,
                                 rng.runif () * 30.0, 
                                 rng.rnorm () * _systemnoise,
                                 this);
            // std::cout << " [" << _state.back ().get_distance () << ", " 
            //     << _state.back ().get_speed () << "," << _state.back ().get_acceleration () << "]";
        }

        double initwt = 1.0 / (double) _N;
        for (auto p = _state.begin (); p != _state.end (); ++p) p->set_weight (initwt);

        _newtrip = false;
        _complete = false;
        bad_sample = false;
        resample = true;

        _current_segment = 0;
        _current_stop = 0;
        _segment_travel_times.clear ();
        _stop_arrival_times.clear ();

        if (_trip == nullptr || _trip->shape () == nullptr) return;
        _segment_travel_times.resize (_trip->shape ()->segments ().size (), 0);
        _stop_arrival_times.resize (_trip->stops ().size (), 0);
    }

    void Vehicle::mutate (RNG& rng)
    {
        if (_newtrip)
        {
            initialize (rng);
            return;
        }

        if (_complete || !valid () || _delta == 0) return;

        // There probably need to be a bunch of checks here ...
        

        // do the transition ("mutation")
        int ncomplete = 0;
        for (auto p = _state.begin (); p != _state.end (); ++p)
        {
            if (p->is_complete ())
            {
                ncomplete++;
            }
            else
            {
                p->travel (_delta, rng);
            }
        }
        if (ncomplete == _state.size ())
        {
            _complete = true;
            return;
        }

#if WRITE_PARTICLES
        std::vector<ShapePt>* path = &(_trip->shape ()->path ());
        for (auto p = _state.begin (); p != _state.end (); ++p)
        {
            p->calculate_likelihood (_position, path, _gpserror);
            std::ostringstream fname;
            fname << "history/vehicle_" << _vehicle_id << "_proposals.csv";
            std::ofstream fout;
            fout.open (fname.str ().c_str (), std::ofstream::app);
            double d (p->get_distance ());
            latlng ppos (_trip->shape ()->coordinates_of (d));
            fout << _timestamp << ","
                << _trip->trip_id () << ","
                << p->get_distance () << ","
                << p->get_speed () << ","
                << p->get_acceleration () << ","
                << std::setprecision(15)
                << p->get_ll () << ","
                << ppos.latitude << "," << ppos.longitude << "\n";
            fout.close ();
        }
#endif

#if SIMULATION
        // PRIOR model eval stuff
        double prior_mse = 0.0; // = SUM [W * d(h(X), Y)^2]
        double dxy;
        latlng hx;
        for (auto p = _state.begin (); p != _state.end (); ++p)
        {
            dxy = p->get_distance ();
            hx = _trip->shape ()->coordinates_of (dxy);
            prior_mse += p->get_weight () * 
                pow(distanceEarth (_position, hx), 2);
        }
#endif
        // update
        select (rng);

        // (re)initialize if the particle sample is bad
        if (bad_sample)
        {
            initialize (rng);
        }
        else
        {
            // NETWORK STUFF
            double dmin = _trip->shape ()->path ().back ().distance;
            for (auto p = _state.begin (); p != _state.end (); ++p)
            {
                if (p->get_distance () < dmin)
                {
                    dmin = p->get_distance ();
                }
            }
            std::vector<ShapeSegment>& segs = _trip->shape ()->segments ();
            int m = find_stop_index (dmin, &(_trip->stops ()));
            int l = find_segment_index (dmin, &segs);


            // update segment travel times for intermediate ones ...
            double tt, ttp, err;
            int n;
            while (_current_segment < m)
            {
                // get the average travel time for particles along that segment
                tt = 0.0;
                n = 0;
                for (auto p = _state.begin (); p != _state.end (); ++p)
                {
                    ttp = p->get_travel_time (_current_segment) * p->get_weight ();
                    if (ttp > 0)
                    {
                        tt += ttp;
                        n++;
                    }
                }
                if (n < _N)
                {
                    _current_segment++;
                    continue;
                }

                err = std::accumulate (_state.begin (), _state.end (), 0.0,
                                       [=](double a, Particle& p) {
                                            return a + p.get_weight () * pow(p.get_travel_time (_current_segment) - tt, 2);
                                       });

                _segment_travel_times.at (_current_segment) = round (tt);
                segs.at (_current_segment).segment->push_data (round (tt), fmax(2.0, err));
                _current_segment++;
            }
            
            // NOTE: need to ignore segment if previous segment travel time is 0
            // (i.e., can't be sure that the current segment travel time is complete)
            
            
        }

#if SIMULATION
        // POSTERIOR model eval stuff
        double posterior_mse = 0.0; // = SUM [W * d(h(X), Y)^2]
        for (auto p = _state.begin (); p != _state.end (); ++p)
        {
            dxy = p->get_distance ();
            hx = _trip->shape ()->coordinates_of (dxy);
            posterior_mse += p->get_weight () * 
                pow(distanceEarth (_position, hx), 2);
        }

        std::ostringstream mename;
        mename << "modeleval/vehicle_" << _vehicle_id << ".csv";
        std::ofstream modeleval;
        modeleval.open (mename.str ().c_str (), std::ofstream::app);
        double atd = _trip->shape ()->distance_of (_position);
        latlng closest_pt = _trip->shape ()->coordinates_of (atd);
        double ctd = distanceEarth (_position, closest_pt);
        modeleval << _vehicle_id 
            << "," << _trip->trip_id ()
            << "," << _timestamp
            << "," << prior_mse 
            << "," << posterior_mse
            << "," << ctd
            << "," << _Neff
            << "," << (resample ? 1 : 0)
            << "\n";
        modeleval.close ();
#endif


#if WRITE_PARTICLES
        for (auto p = _state.begin (); p != _state.end (); ++p)
        {
            p->calculate_likelihood (_position, path, _gpserror);
            std::ostringstream fname;
            fname << "history/vehicle_" << _vehicle_id << ".csv";
            std::ofstream fout;
            fout.open (fname.str ().c_str (), std::ofstream::app);
            double d (p->get_distance ());
            latlng ppos (_trip->shape ()->coordinates_of (d));
            fout << _timestamp << ","
                << _trip->trip_id () << ","
                << std::setprecision(15)
                << _position.latitude << "," << _position.longitude << ","
                << std::setprecision(6)
                << p->get_distance () << ","
                << p->get_speed () << ","
                << p->get_acceleration () << ","
                << std::setprecision(15)
                << p->get_ll () << ","
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
        }
#endif


    }

    void Vehicle::select (RNG& rng)
    {
        bad_sample = true;
        resample = false;

        // calculate loglikelihood of particles
        double sumlh = 0.0;
        // threshold of 100m
        double threshold = log (0.5) - 0.5 * exp (2.0 * (log (100.0) - log (_gpserror)));
        double maxlh = threshold;
        std::vector<ShapePt>* path = &(_trip->shape ()->path ());
        double plh;
        double sumwt = 0.0;
        for (auto p = _state.begin (); p != _state.end (); ++p)
        {
            p->calculate_likelihood (_position, path, _gpserror);
            plh = p->get_ll ();
            sumlh += exp (plh);
            if (plh > maxlh) maxlh = plh;

            p->set_weight (p->get_weight () * exp (plh));
            sumwt += p->get_weight ();
        }

        // if (maxlh < threshold) return;
        
        if (sumwt == 0) return;

        // normalize weights
        std::vector<double> wt;
        wt.reserve (_N);
        for (auto p = _state.begin (); p != _state.end (); ++p)
        {
            p->set_weight (p->get_weight () / sumwt);
            wt.push_back (p->get_weight ());
        }


        // some condition for resampling (1/wt^2 < N? or something ...)
        sumwt = std::accumulate (wt.begin (), wt.end (), 0.0);
        if (sumwt < 0.9) return;

        {
            double sumwt2 = std::accumulate(wt.begin (), wt.end (), 0.0,
                                            [](double a, double b) {
                                                return a + pow(b, 2);
                                            });
            _Neff = 1.0 / sumwt2;
            // std::cout << "\n sumwt2 = " << sumwt2 <<", Neff = " << _Neff;
        }

        bad_sample = false;
        
        if (_Neff >= 1000) return;
        resample = true;

        // resample u \in [0, 1)
        double u;
        unsigned int j;
        std::vector<Particle> _newstate;
        _newstate.reserve (_N);
        for (int i=0; i<_N; ++i)
        {
            u = rng.runif () * sumwt;
            for (j=0; j<_N; ++j)
            {
                if (u < wt[j])
                {
                    _newstate.emplace_back (_state[j]);
                    break;
                }
                // if not less, subtract wt from u
                u -= wt[j];
            }
        }

        // std::cout << "\n + vehicle " << _vehicle_id << ":\n";
        // for (auto p = _state.begin (); p != _state.end (); ++p)
        //     std::cout << "[" << p->get_distance () << ", " << p->get_speed () << "] ; ";
        _state.clear ();
        _state = std::move(_newstate);
        // std::cout << "\n";
        // for (auto p = _state.begin (); p != _state.end (); ++p)
        //     std::cout << "[" << p->get_distance () << ", " << p->get_speed () << "] ; ";

        double initwt = pow(_N, -1);
        for (auto p = _state.begin (); p != _state.end (); ++p) p->set_weight (initwt);
        
    }

    void Vehicle::predict_etas (RNG& rng)
    {
        if (!valid ()) return;

#if VERBOSE == 2
        Timer timer;
        std::cout << "- vehicle " << _vehicle_id << " - predicting etas";
#endif
        for (auto p = _state.begin (); p != _state.end (); ++p) 
        {
            p->predict_etas (rng);
            // std::cout << "\n";
            // for (auto eta : p->get_arrival_times ()) std::cout << eta << ", ";
            // std::cout << "\n VERSUS ... ";
            // for (int i=0; i<p->get_arrival_times ().size (); ++i) 
            //     std::cout << p->get_arrival_time (i) << ",";
        }
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
        if (!valid ()) return etas;
        // std::cout << "\n Vehicle " << _vehicle_id << " =============================";
        std::vector<uint64_t> etam;
        std::vector<double> wts;
        etam.reserve (_state.size ());
        wts.reserve (_state.size ());
        double sumwt;
        for (int i=0; i<M; ++i)
        {
            // need to center each particle's arrival time
            // std::cout << "\n   [" << i << "]: ";
            int tarr = 0;
            int ni = 0;
            etam.clear ();
            wts.clear ();

            for (auto p = _state.begin (); p != _state.end (); ++p)
            {
                // std::cout << p->get_arrival_time (i) << ", ";
                if (p->get_arrival_time (i) > 0)
                {
                    tarr += p->get_weight () * (p->get_arrival_time (i) - _timestamp);
                    etam.push_back (p->get_arrival_time (i));
                    wts.push_back (p->get_weight ());
                }
            }
            if (wts.size () != _N) continue;
            sumwt = std::accumulate (wts.begin (), wts.end (), 0.0);
            if (sumwt != 1.0) continue;

            etas.at (i).stop_id = stops.at (i).stop->stop_id ();
            etas.at (i).estimate = _timestamp + tarr;
            // generate quantiles [0, 5, 50, 95, 100] -> [0, 50, 500, 950, 999]
            std::sort (etam.begin (), etam.end ());
            etas.at (i).quantiles.emplace_back (0.0, etam.front ());
            int qi = 0;
            double cwt = 0.0;
            while (cwt <= 0.05 && qi < wts.size ())
            {
                cwt += wts.at (qi);
                qi++;
                // std::cout << ".";
            }
            etas.at (i).quantiles.emplace_back (5.0, etam.at (qi));
            while (cwt <= 0.5 && qi < wts.size ())
            {
                cwt += wts.at (qi);
                qi++;
                // std::cout << ".";
            }
            etas.at (i).quantiles.emplace_back (50.0, etam.at (qi));
            while (cwt <= 0.95 && qi < wts.size ())
            {
                cwt += wts.at (qi);
                qi++;
                // std::cout << ".";
            }
            etas.at (i).quantiles.emplace_back (95.0, etam.at (qi));
            etas.at (i).quantiles.emplace_back (100.0, etam.back ());
        }
        return etas;
    }

    double Vehicle::distance ()
    {
        double distance = 0.0;
        for (auto p = _state.begin (); p != _state.end (); ++p) 
            distance += p->get_weight () * p->get_distance ();
        return distance;
    }

    double Vehicle::speed ()
    {
        double speed = 0.0;
        for (auto p = _state.begin (); p != _state.end (); ++p) 
            speed += p->get_weight () * p->get_speed ();
        return speed;
    }

    int Vehicle::progress ()
    {
        double d = distance ();
        double dmax = _trip->shape ()->path ().back ().distance;
        return (100 * d / dmax) + 0.5;
    }



    void Particle::travel (unsigned delta, RNG& rng)
    {
        if (!vehicle || !vehicle->trip () || !vehicle->trip ()->shape ()) return;
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
        stops = &(vehicle->trip ()->stops ());
        int M (stops->size ());
        unsigned int m (find_stop_index (distance, stops));
        if (m == M-1) 
        {
            distance = Dmax;
            complete = true;
            return;
        }

        double next_stop_d = stops->at (m+1).distance;
        
        // get SEGMENTS
        std::vector<ShapeSegment>* segments;
        segments = &(vehicle->trip ()->shape ()->segments ());
        int L (segments->size ());
        unsigned int l (find_segment_index (distance, segments));
        double next_segment_d;
        next_segment_d = (l+1 >= L-1) ? Dmax : segments->at (l+1).distance;
        
        // allow vehicle to remain stationary if at a stop:
        if (distance == stops->at (m).distance &&
            rng.runif () < 0.05)
        {
            double w = - log (rng.runif ()) * delta;
            delta = fmax (0, delta - round (w));
            // we don't want this to affect the speed
        }
        else if (distance == segments->at (l).distance &&
                 rng.runif () < 0.05)
        {
            double w = - log (rng.runif ()) * delta;
            delta = fmax (0, delta - round (w));
            if (tt.at (l) >= 0)
                tt.at (l) = tt.at (l) + 1;
        }
        else if (rng.runif () < 0.05)
        {
            // a very small chance for particles to remain stationary
            speed = 0.0;
            // when /not/ at a bus stop, set speed to 0 and wait
            double w = - log (rng.runif ()) * delta;
            delta = fmax (0.0, delta - round (w));
            if (tt.at (l) >= 0)
                tt.at (l) = tt.at (l) + 1;
            // then the bus needs to accelerate back up to speed ... for how many seconds?
            accelerating = 5.0 + rng.runif () * 10.0;
            acceleration = 2.0 + rng.rnorm () * vehicle->system_noise ();
        }


        double speed_mean = 10.0;
        double speed_sd = 100.0;
        if (segments->at (l).segment->travel_time () > 0 &&
            segments->at (l).segment->uncertainty () > 0) 
        {
            speed_mean = segments->at (l).segment->get_speed ();
            speed_sd = - segments->at (l).segment->length () / pow (speed_mean, 2) *
                segments->at (l).segment->uncertainty ();
            speed_sd = pow(speed_sd, 0.5);
        }
        while (distance < Dmax && delta > 0.0)
        {
            // add system noise to acceleration to ensure speed remains in [0, 30]
            double accel_prop (-100.0);
            double n = 0;
            while (speed + accel_prop < 0.0 || speed + accel_prop > 30.0)
            {
                accel_prop = rng.rnorm () * vehicle->system_noise () * 
                    (1.0 + (double)n / 100.0);
                n++;
                if (accelerating > 0.0)
                {
                    accel_prop += acceleration;
                    accelerating--;
                }
            }

            double v = fmax (0, fmin (30, speed + acceleration));
            double vstar = speed + accel_prop;
            double alpha = (pow (v - speed_mean, 2) - pow(vstar - speed_mean, 2)) / (2 * pow (speed_sd, 2));
            alpha = fmin (0, alpha);
            if (rng.runif () < exp (alpha))
            {
                acceleration = accel_prop;
                speed = vstar;
            }
            else
            {
                speed = v;
            }

            distance += speed;
            delta--;
            if (tt.at (l) >= 0)
                tt.at (l) = tt.at (l) + 1;
            
            if (l < L-1 && distance >= next_segment_d)
            {
                // reaching intersection ... 
                l++;
                tt.at (l) = 0;
                next_segment_d = (l+1 >= L-1) ? Dmax : segments->at (l+1).distance;
                if (segments->at (l).segment->travel_time () > 0 &&
                    segments->at (l).segment->uncertainty () > 0) 
                {
                    speed_mean = segments->at (l).segment->get_speed ();
                    speed_sd = - segments->at (l).segment->length () / pow (speed_mean, 2) *
                        segments->at (l).segment->uncertainty ();
                    speed_sd = pow(speed_sd, 0.5);
                }
                else
                {
                    speed_mean = 10.0;
                    speed_sd = 100.0;
                }
            }

            if (distance >= next_stop_d)
            {
                // about to reach a stop ... slow? stop? just drive past?
                m++; // the stop we are about to reach
                at.at (m) = vehicle->timestamp () - delta;
                if (m == M-1) 
                {
                    distance = next_stop_d;
                    break;
                }
                if (rng.runif () < vehicle->pr_stop ())
                {
                    // stop dwell time ~ Exp(tau = 10)
                    double dwell = vehicle->gamma () - vehicle->dwell_time () * log (rng.runif ());
                    delta = fmax(0, delta - dwell);
                    distance = next_stop_d;
                    speed = 0.0;
                    accelerating = 5.0 + rng.runif () * 10.0;
                    acceleration = 2.0 + rng.rnorm () * vehicle->system_noise ();
                }
                next_stop_d = stops->at (m+1).distance;
                continue;
            }
        }
    }

    void Particle::predict_etas (RNG& rng)
    {
        // std::cout << "\n > ";
        if (complete) return;
        // std::cout << "| ";

        // get STOPS
        std::vector<StopTime>* stops;
        if (!vehicle || !vehicle->trip ()) return;
        stops = &(vehicle->trip ()->stops ());
        int M (stops->size ());
        if (M == 0) return;
        at.clear ();
        at.resize (M, 0);
        unsigned int m (find_stop_index (distance, stops));
        // std::cout << " @" << m << " > ";
        if (m == M-1) 
        {
            return;
        }
        double Dmax = stops->back ().distance;
        
        if (!vehicle->trip ()->shape ()) return;
        std::vector<ShapeSegment>* segments;
        segments = &(vehicle->trip ()->shape ()->segments ());
        int L (segments->size ());
        unsigned int l (find_segment_index (distance, segments));
        double next_segment_d;
        next_segment_d = (l+1 >= L-1) ? Dmax : segments->at (l+1).distance;

        double dcur = distance;
        double dnext;
        uint64_t t0 = vehicle->timestamp ();
        // std::cout << t0 << " > ";
        int etat;
        double vel;
        vel = segments->at (l).segment->sample_speed (rng);
        if (vel == 0.0 && speed >= 0.0) vel = speed;
        while (vel <= 0 || vel > 30)
        {
            vel = rng.rnorm () * 8.0 + 15.0;
        }
        while (m < M-1)
        {
            m++; // `next` stop index
            dnext = stops->at (m).distance;
            etat = 0;
            // std::cout << " [";
            while (next_segment_d < dnext && l < L-1)
            {
                // time to get to end of segment
                etat += (next_segment_d - dcur) / vel;
                dcur = next_segment_d;
                l++;
                next_segment_d = (l+1 >= L-1) ? Dmax : segments->at (l+1).distance;
                vel = segments->at (l).segment->sample_speed (rng);
                if (vel == 0.0 && speed >= 0.0) vel = speed;
                while (vel <= 0.0 || vel > 30)
                {
                    vel = rng.rnorm () * 8.0 + 15.0;
                }
                // std::cout << vel << "; ";
            }
            // std::cout << "] ";
            etat += (dnext - dcur) / vel;
            at.at (m) = t0 + etat; // makes no sense because speeds are noise
            dcur = dnext;
            // and add some dwell time
            t0 = at.at (m);
            if (rng.runif () < vehicle->dwell_time ())
            {
                t0 += vehicle->gamma () - vehicle->dwell_time () * log (rng.runif ());
            }
            // std::cout << "(" << m << ") " << at.at (m) << ", ";
        }
    }

    void
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
    }

    void Particle::set_weight (double w)
    {
        weight = w;
    }

}; // namespace Gtfs
