#include "paa.h"


#include <string>
#include <algorithm>
using std::vector;

float paa::get_mean(const float* const first, const float* const last)
{
  if (first > last) return 0.0;
  float sum = 0; 

  for (int i=0; i<=last-first; ++i) {
    sum += *(first+i);
  }
  return sum / (last - first + 1);
}

vector<float> paa::paa(const vector<float>& series, unsigned int interval_size)
{
  if (series.size() == 0) return {};

  vector<float> paa_vec;

  int num_full_intervals = (series.size() / interval_size);
  for (int i=0; i<num_full_intervals; ++i) {
    int start_pos = i*num_full_intervals;
    int end_pos = start_pos + interval_size;
    paa_vec.emplace_back( get_mean(series.data() + start_pos, series.data() + end_pos) );
  }
  
  if (series.size() % interval_size != 0) {
    int start_pos = num_full_intervals * interval_size;
    int end_pos = start_pos + (series.size() % interval_size) - 1;
    paa_vec.emplace_back( get_mean( series.data() + start_pos, series.data() + end_pos) );
  }
  return paa_vec;
}
