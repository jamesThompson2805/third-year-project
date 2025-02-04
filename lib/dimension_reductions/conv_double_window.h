#ifndef CONV_DOUBLE_WINDOW
#define CONV_DOUBLE_WINDOW

#include <vector>
#include <tuple>

#include "pla.h"

namespace c_d_w {
  std::vector<std::tuple<DoublePair,unsigned int>> conv_pla(const std::vector<double>& s, unsigned int num_params, const std::vector<double>& l, const std::vector<double>& r);
}

#endif
