#ifndef DOUBLE_WINDOW_H
#define DOUBLE_WINDOW_H

#include <vector>
#include <tuple>

#include "pla.h"

namespace d_w {
  std::vector<std::tuple<double,unsigned int>> simple_paa(const std::vector<double>& s, unsigned int num_params, unsigned int lw_size, unsigned int rw_size);
  std::vector<std::tuple<double,unsigned int>> y_proj_paa(const std::vector<double>& s, unsigned int num_params, unsigned int lw_size, unsigned int rw_size);

  std::vector<std::tuple<DoublePair,unsigned int>> simple_pla(const std::vector<double>& s, unsigned int num_params, unsigned int lw_size, unsigned int rw_size);
  std::vector<std::tuple<DoublePair,unsigned int>> y_proj_pla(const std::vector<double>& s, unsigned int num_params, unsigned int lw_size, unsigned int rw_size);
}

#endif
