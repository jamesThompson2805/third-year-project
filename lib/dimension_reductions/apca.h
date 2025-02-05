#include <vector>
#include <tuple>

namespace apca {
  std::vector<std::tuple<double, unsigned int>> apca(const std::vector<double>& s, unsigned int num_params);
}
