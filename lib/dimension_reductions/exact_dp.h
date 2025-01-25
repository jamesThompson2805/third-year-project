#ifndef EXACT_H
#define EXACT_H

#include <vector>
#include <tuple>

#include "pla.h"

namespace exact_dp { 
std::vector< std::tuple< double, unsigned int> > min_mse_paa( const std::vector<double>&, unsigned int num_params);
std::vector< std::tuple< DoublePair, unsigned int> > min_mse_pla( const std::vector<double>&, unsigned int num_params);

std::vector< std::tuple< double, unsigned int> > min_maxdev_paa( const std::vector<double>&, unsigned int num_params);
std::vector< std::tuple< DoublePair, unsigned int> > min_maxdev_pla( const std::vector<double>&, unsigned int num_params);
};

#endif

