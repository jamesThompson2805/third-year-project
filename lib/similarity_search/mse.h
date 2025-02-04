#ifndef MSE_H
#define MSE_H

#include <vector>

namespace mse {
  double se_between_seq(const std::vector<double>& s1, const std::vector<double>& s2);
  double mse_between_seq(const std::vector<double>& s1, const std::vector<double>& s2);
  double maxdev_between_seq(const std::vector<double>& s1, const std::vector<double>& s2);
};
#endif
