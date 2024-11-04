#ifndef DAC_CURVE_FITTING_H
#define DAC_CURVE_FITTING_H

#include <vector>
#include <tuple>

#include "pla.h"

namespace dac_curve_fitting {

  std::vector<std::tuple<double, unsigned int>> dac_constant( const std::vector<double>& s, double epsilon); 
  std::vector<std::tuple<double, unsigned int>> dac_constant_intervals( const std::vector<double>& s, unsigned int num_intervals); 

  std::vector<std::tuple<DoublePair, unsigned int>> dac_linear( const std::vector<double>& s, double epsilon);

}

#endif
