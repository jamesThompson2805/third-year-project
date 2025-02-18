#include "error_measures.h"

#include <algorithm>
#include <cmath>
using std::vector;

double error_measures::se_between_seq(const vector<double>& s1, const vector<double>& s2)
{
  double se = 0;
  for (int i=0; i<std::min( s1.size(), s2.size() ); ++i) {
    se += (s1[i] - s2[i]) * (s1[i] - s2[i]);
  }
  return se;
}
double error_measures::se_between_ptrs(const double* const s1_start, const double* const s1_end, const double* const s2_start, const double* const s2_end)
{
  double se = 0;
  for (int i=0; i<=std::min( s1_end-s1_start, s2_end-s2_start ); ++i) {
    se += ( s1_start[i] - s2_start[i]) * (s1_start[i] - s2_start[i]);
  }
  return se;
}
double error_measures::mse_between_seq(const vector<double>& s1, const vector<double>& s2)
{
  return error_measures::se_between_seq(s1, s2) / std::min( s1.size(), s2.size()) ;
}
double error_measures::l2_between_seq(const vector<double>& s1, const vector<double>& s2)
{
  return std::sqrt(error_measures::se_between_seq(s1, s2));
}

double error_measures::maxdev_between_seq(const vector<double>& s1, const vector<double>& s2)
{
  double max_dev = -1.0;
  for (int i=0; i<std::min( s1.size(), s2.size() ); ++i) {
    if ( double dev = std::abs( s1[i] - s2[i] ); dev > max_dev)
      max_dev = dev;
  }
  return max_dev;
}
