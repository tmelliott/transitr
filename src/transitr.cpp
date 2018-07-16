#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector test ()
{
  NumericVector x (1);
  x[0] = 5;

  return x;
}
