#ifndef EXACT_H
#define EXACT_H

/**
 * @file exact_dp.h contains methods for computing the optimal partitions for constant and linear approximations
 */

#include <vector>
#include <tuple>

#include "pla.h"

/**
 * @brief exact_dp namespace contains all methods of exact_dp.h
 */
namespace exact_dp { 
  /**
   * @brief min_l2_paa finds optimal partition for paa under euclidean distance
   * @param s is series to compress
   * @param num_params is the target dimension
   * @return sorted array of segments 
   */
std::vector< std::tuple< double, unsigned int> > min_l2_paa( const Seqd& s, unsigned int num_params);
  /**
   * @brief min_l2_pla finds optimal partition for pla under euclidean distance
   * @param s is series to compress
   * @param num_params is the target dimension
   * @return sorted array of segments 
   */
Seqddt min_l2_pla( const Seqd&, unsigned int num_params);

  /**
   * @brief min_maxdev_paa finds optimal partition for paa under maximum deviation
   * @param s is series to compress
   * @param num_params is the target dimension
   * @return sorted array of segments 
   */
std::vector< std::tuple< double, unsigned int> > min_maxdev_paa( const Seqd& s, unsigned int num_params);
  /**
   * @brief min_maxdev_pla finds optimal partition for pla under maximum deviation
   * @param s is series to compress
   * @param num_params is the target dimension
   * @return sorted array of segments 
   */
Seqddt min_maxdev_pla( const Seqd&, unsigned int num_params);

};

#endif

