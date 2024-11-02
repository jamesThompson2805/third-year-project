#ifndef MSE_H
#define MSE_H

#include <vector>

namespace mse {
  float mse_between_seq(const std::vector<float>& s1, const std::vector<float>& s2);
};
#endif
