#include "paa.h"

using std::vector;

double paa::get_mean(const double* const first, const double* const last)
{
  if (first > last) return 0.0;
  double sum = 0; 

  for (int i=0; i<=last-first; ++i) {
    sum += *(first+i);
  }
  return sum / (last - first + 1);
}

vector<double> paa::paa(const vector<double>& series, unsigned int interval_size)
{
  if (series.size() == 0) return {};

  vector<double> paa_vec;

  int num_full_intervals = (series.size() / interval_size);
  for (int i=0; i<num_full_intervals; ++i) {
    int start_pos = i*interval_size;
    int end_pos = start_pos + interval_size - 1;
    paa_vec.emplace_back( get_mean(series.data() + start_pos, series.data() + end_pos) );
  }
  
  if (series.size() % interval_size != 0) {
    int start_pos = num_full_intervals * interval_size;
    int end_pos = start_pos + (series.size() % interval_size) - 1;
    paa_vec.emplace_back( get_mean( series.data() + start_pos, series.data() + end_pos) );
  }
  return paa_vec;
}

double paa::paa_mse(const vector<double>& series, unsigned int interval_size)
{
  vector<double> paa_series = paa::paa(series, interval_size);
  double mse = 0;
  for (int i=0; i<series.size(); ++i) {
    double diff = (series[i] - paa_series[ i / interval_size ]);
    mse += diff * diff;
  }
  return mse / series.size();
}
