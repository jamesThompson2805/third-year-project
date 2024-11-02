#ifndef PLA_H
#define PLA_H

#include <array>
#include <vector>

typedef std::array<float, 2> FloatPair;

namespace pla {
std::vector<std::vector<float>> chunk_series(std::vector<float> series, unsigned int chunk_size);


FloatPair regression(const float* const start, const float* const end);

// w is num items compressed to a linear function
std::vector<FloatPair> sliding_window_regression( const std::vector<float>& series, unsigned int w);

std::vector<FloatPair> chunk_regression( const std::vector<float>& series, unsigned int interval_size);

float pla_mse(const std::vector<float>& series, unsigned int interval_size);

};

#endif
