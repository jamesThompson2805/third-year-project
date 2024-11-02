#ifndef DAC_CURVE_FITTING_H
#define DAC_CURVE_FITTING_H

#include <vector>
#include <tuple>

namespace dac_curve_fitting {

  std::vector<std::tuple<float, unsigned int>> dac_constant( const std::vector<float>& s, float epsilon); 

  std::vector<std::tuple<float, float, unsigned int>> dac_linear( const float* const s_start, const float* const s_end, unsigned int max_num_intervals);

}

#endif
