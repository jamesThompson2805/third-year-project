#ifndef PAA_H
#define PAA_H
#include <vector>

namespace paa {
  double get_mean(const double* const first, const double* const last);

  std::vector<double> paa(const std::vector<double>& series, unsigned int num_params);

  double paa_mse(const std::vector<double>& series, unsigned int num_params); 
}

#endif
