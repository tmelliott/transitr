#ifndef DISTS_H
#define DISTS_H

#include <math.h>

// HEADER ONLY

namespace stats {
  double normal_pdf_log (double x, double m, double s)
  {
    // return - 2 * log (2 * M_PI * pow (s, 2.)) - pow (x - s, 2.) / 2 / pow (s, 2.);
    return - 2 * log (2) - 2 * log (M_PI) - 4 * log (s) -
      pow (x - s, 2.) * 0.5 * pow (s, -2.);
  }
}

#endif