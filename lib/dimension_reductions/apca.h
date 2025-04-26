#ifndef APCA_H
#define APCA_H
#include <vector>
#include <tuple>

/**
 * @file apca.h is file containing the APCA dimension reduction technique
 */

/**
 * @brief apca namespace contains only the apca approximation
 */
namespace apca {
  /**
   * @brief apca function converts a time series to an APCA form
   * @param s is a series to convert
   * @param num_params is the dimension the approximation will occupy
   * @return an array of pairs of a value and index where the value is the mean on the segment ending at the index, the array is sorted
   */
  std::vector<std::tuple<double, unsigned int>> apca(const std::vector<double>& s, unsigned int num_params);
}

#endif
