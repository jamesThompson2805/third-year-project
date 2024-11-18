#ifndef DAC_CURVE_FITTING_H
#define DAC_CURVE_FITTING_H

//Algorithms used and adapted from paper "Approximate Queries and Representations for Large Data Sequences"
//
#include <vector>
#include <tuple>

#include "pla.h"

namespace dac_curve_fitting {

  std::vector<std::tuple<double, unsigned int>> dac_constant( const std::vector<double>& s, double epsilon); 
  std::vector<std::tuple<double, unsigned int>> dac_constant_params( const std::vector<double>& s, unsigned int num_params); 

  std::vector<std::tuple<DoublePair, unsigned int>> dac_linear( const std::vector<double>& s, double epsilon);
  std::vector<std::tuple<double, unsigned int>> dac_linear_params( const std::vector<double>& s, unsigned int num_params); 

}

#endif
