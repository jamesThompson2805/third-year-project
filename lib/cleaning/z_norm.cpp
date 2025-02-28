#include "z_norm.h"

#include <algorithm>
#include <numeric>
#include <cmath>

#include <iostream>

using std::vector;

void z_norm::z_normalise(vector<double> &series)
{
  double n = series.size();
  double mean = std::accumulate(series.begin(), series.end(), 0.0) / n;  
  double var = 0;
  for (const double &f : series) {
    var += (f-mean)*(f-mean); 
  }
  var /= n;
  double stddev = std::sqrt( var );

  std::transform(series.begin(), series.end(), series.begin(), [&](double& f){ return (f - mean) / stddev; } );
}
