#ifndef PAA_H
#define PAA_H
#include <vector>
#include <tuple>

namespace paa {
  double get_mean(const double* const first, const double* const last);

  std::vector<double> paa(const std::vector<double>& series, unsigned int num_params);

  double paa_mse(const std::vector<double>& series, unsigned int num_params); 

  std::vector<double> paa_to_seq(const std::vector<double>, unsigned int int_size);
  std::vector<double> apca_to_seq(const std::vector< std::tuple<double, unsigned int>>);
}

#endif
