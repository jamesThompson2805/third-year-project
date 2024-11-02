#include "mse.h"

#include <algorithm>
using std::vector;

float mse::mse_between_seq(const vector<float>& s1, const vector<float>& s2)
{
  float mse = 0;
  for (int i=0; i<std::min( s1.size(), s2.size() ); ++i) {
    mse += (s1[i] - s2[i]) * (s1[i] - s2[i]);
  }
  return mse;
}
