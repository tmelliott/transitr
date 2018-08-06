#include "rng.h"

RNG::RNG () : uniform (0.0, 1.0), normal (0.0, 1.0)
{
}

RNG::RNG (unsigned int seed) : uniform (0.0, 1.0), normal (0.0, 1.0)
{
    set_seed (seed);
}

void RNG::set_seed (unsigned int seed)
{
    gen.seed (seed);
}

double RNG::runif ()
{
    return uniform (gen);
}

double RNG::rnorm ()
{
    return normal (gen);
}
