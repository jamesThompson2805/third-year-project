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

vector<double> paa::paa(const vector<double>& series, unsigned int num_params)
{
  if (series.size() == 0) return {};

  vector<double> paa_vec;
  unsigned int interval_size = (series.size() / num_params) + (series.size() % num_params != 0);

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

double paa::paa_mse(const vector<double>& series, unsigned int num_params)
{
  vector<double> paa_series = paa::paa(series, num_params);
  unsigned int interval_size = (series.size() / num_params) + (series.size() % num_params != 0);
  double mse = 0;
  for (int i=0; i<series.size(); ++i) {
    double diff = (series[i] - paa_series[ i / interval_size ]);
    mse += diff * diff;
  }
  return mse / series.size();
}

vector<double> paa::paa_to_seq(const std::vector<double> paa_s, unsigned int int_size)
{
  vector<double> seq;
  for (const auto& d : paa_s) {
    for (int i=0; i<int_size; ++i)
      seq.emplace_back(d);
  }
  return seq;
}

vector<double> paa::apca_to_seq(const std::vector<std::tuple<double, unsigned int>> apca_seq)
{
  vector<double> seq;
  int curr_index = 0;
  for (const auto& [d, i] : apca_seq) {
    for (; curr_index <= i; ++curr_index)
      seq.emplace_back(d);
  }
  return seq;
}
