#ifndef PLA_H
#define PLA_H

#include <array>
#include <vector>
#include <tuple>

typedef std::array<double, 2> DoublePair;

namespace pla {
std::vector<std::vector<double>> chunk_series(std::vector<double> series, unsigned int chunk_size);


DoublePair regression(const double* const start, const double* const end);

// w is num items compressed to a linear function
std::vector<DoublePair> sliding_window_regression( const std::vector<double>& series, unsigned int w);

std::vector<DoublePair> chunk_regression( const std::vector<double>& series, unsigned int num_params);

double pla_mse(const std::vector<double>& series, unsigned int interval_size);

std::vector<double> pla_to_seq(const std::vector<DoublePair>, unsigned int int_size);
std::vector<double> apla_to_seq(const std::vector< std::tuple<DoublePair, unsigned int>>);

};

#endif
