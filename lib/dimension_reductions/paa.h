#include <vector>

namespace paa {
  float get_mean(const float* const first, const float* const last);

  std::vector<float> paa(const std::vector<float>& series, unsigned int interval_size);
}
