#include "mse.h"

#include <algorithm>
#include <cmath>
using std::vector;

double mse::se_between_seq(const vector<double>& s1, const vector<double>& s2)
{
  double se = 0;
  for (int i=0; i<std::min( s1.size(), s2.size() ); ++i) {
    se += (s1[i] - s2[i]) * (s1[i] - s2[i]);
  }
  return se;
}
double mse::mse_between_seq(const vector<double>& s1, const vector<double>& s2)
{
  return std::sqrt( mse::se_between_seq(s1, s2) ) / std::min( s1.size(), s2.size()) ;
}

double mse::maxdev_between_seq(const vector<double>& s1, const vector<double>& s2)
{
  double max_dev = -1.0;
  for (int i=0; i<std::min( s1.size(), s2.size() ); ++i) {
    if ( double dev = std::abs( s1[i] - s2[i] ); dev > max_dev)
      max_dev = dev;
  }
  return max_dev;
}
