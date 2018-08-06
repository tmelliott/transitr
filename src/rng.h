#ifndef RNG_H
#define RNG_H

#include <random>

class RNG
{
private:
    std::mt19937_64 gen;

    std::uniform_real_distribution<double> uniform;
    std::normal_distribution<double> normal;


public:
    RNG ();
    RNG (unsigned int seed);

    void set_seed (unsigned int seed);

    double runif ();
    double rnorm ();

};


#endif
