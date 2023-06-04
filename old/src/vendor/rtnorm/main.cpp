//  Example for using rtnorm
//  
//  Copyright (C) 2012 Guillaume Dollé, Vincent Mazet (LSIIT, CNRS/Université de Strasbourg)
//  Licence: GNU General Public License Version 2
//  see http://www.gnu.org/licenses/old-licenses/gpl-2.0.txt
//
//  Depends: LibGSL
//  OS: Unix based system


#include <iostream>
#include <gsl/gsl_rng.h>

#include "rtnorm.hpp"


int main()
{
  double a = 1;                 // Left bound
  double b = 9;                 // Right bound
  double mu = 2;                // Mean
  double sigma = 3;             // Standard deviation
  std::pair<double, double> s;  // Output argument of rtnorm
  int K = 1e5;                  // Number of random variables to generate

  //--- GSL random init ---
  gsl_rng_env_setup();                          // Read variable environnement
  const gsl_rng_type* type = gsl_rng_default;   // Default algorithm 'twister'
  gsl_rng *gen = gsl_rng_alloc (type);          // Rand generator allocation

  //--- generate and display the random numbers ---
  //std::cout<<"# x p(x)"<<std::endl;
    std::cout<<a<<" "<<b<<std::endl;
    std::cout<<mu<<" "<<sigma<<std::endl;
  for(int k=0; k<K; k++)
  {
    s = rtnorm(gen,a,b,mu,sigma);
    std::cout<<s.first<<" "<<s.second<<std::endl;
  }

  gsl_rng_free(gen);                            // GSL rand generator deallocation

  return 0;
}

