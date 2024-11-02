#include "z_norm.h"

#include <algorithm>
#include <numeric>
#include <cmath>

#include <iostream>

using std::vector;

void z_norm::z_normalise(vector<float> &series)
{
  float n = series.size();
  float mean = std::accumulate(series.begin(), series.end(), 0.0) / n;  
  float var = 0;
  for (const float &f : series) {
    var += (f-mean)*(f-mean); 
  }
  float stddev = std::sqrt( var );

  std::cout << mean << " : " << stddev << std::endl;

  std::transform(series.begin(), series.end(), series.begin(), [&](float& f){ return (f - mean) / stddev; } );
}
