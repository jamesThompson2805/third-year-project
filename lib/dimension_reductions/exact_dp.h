#ifndef EXACT_H
#define EXACT_H

#include <vector>
#include <tuple>

#include "pla.h"

namespace exact_dp { 
std::vector< std::tuple< double, unsigned int> > exact_paa( const std::vector<double>&, unsigned int num_params);
std::vector< std::tuple< DoublePair, unsigned int> > exact_pla( const std::vector<double>&, unsigned int num_params);
};

#endif

