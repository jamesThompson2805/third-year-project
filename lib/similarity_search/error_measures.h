#ifndef MSE_H
#define MSE_H

#include <vector>

namespace error_measures {
  double se_between_seq(const std::vector<double>& s1, const std::vector<double>& s2);
  double mse_between_seq(const std::vector<double>& s1, const std::vector<double>& s2);
  double l2_between_seq(const std::vector<double>& s1, const std::vector<double>& s2);
  double maxdev_between_seq(const std::vector<double>& s1, const std::vector<double>& s2);
};
#endif
