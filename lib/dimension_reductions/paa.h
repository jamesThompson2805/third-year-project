#ifndef PAA_H
#define PAA_H

/**
 * @file paa.h contains the paa approximation
 */

#include <vector>
#include <tuple>

/**
 * @brief paa namespace contains methods to convert sequence to paa, convert paa to sequence and convert adaptive paa to sequence
 */
namespace paa {
  /**
   * @brief get_mean is fast method to find mean of an array, passed by a start and end pointer
   * @param first is pointer to first element
   * @param last is pointer to last element
   * @return mean of all elements
   */
  double get_mean(const double* const first, const double* const last);

  /**
   * @brief paa calculates the paa of a series
   * @param series is the series to compress
   * @param num_params is the target dimension
   * @return sorted array representing the mean of each segment
   */
  std::vector<double> paa(const std::vector<double>& series, unsigned int num_params);

  /**
   * @brief paa_mse calculates paa of series and finds the mean squared error of the approximation
   * @param series is the series to compress and analyse compression against
   * @param num_params is target dimension of the compression
   * @return the mean squared error of approximation to original
   */
  double paa_mse(const std::vector<double>& series, unsigned int num_params); 

  /**
   * @brief paa_to_seq takes a paa approximation and its segment size and returns the full size approximation
   * @param paa_s is the compressed series
   * @param int_size is the size of a segment
   * @return a series representing the uncompressed paa_s
   */
  std::vector<double> paa_to_seq(const std::vector<double> paa_s, unsigned int int_size);
  /**
   * @brief apca_to_seq takes a apca approximation and returns the uncompressed approximation
   * @param apca_s is the compressed series
   * @return a series representing the uncompressed apca_s
   */
  std::vector<double> apca_to_seq(const std::vector< std::tuple<double, unsigned int>> apca_s);
}

#endif
