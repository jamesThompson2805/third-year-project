#ifndef DOUBLE_WINDOW_H
#define DOUBLE_WINDOW_H

#include "pla.h"

/**
 * @file double_window.h is header file implementing the double sliding window algorithms
 * These implementations for Average Value and Interval Projection are unoptimised
 */

/**
 * @brief d_w namespace holds unoptimised double window DRTs
 */
namespace d_w {
  /**
   * @brief simple_pla is the average value approximation
   * @param s is the series to compress
   * @param num_params is the target dimension for the compression
   * @param lw_size is the size of the first window
   * @param rw_size is the size of the second window
   * @return the Adaptive PLA representation
   */
  Seqddt simple_pla(const Seqd& s, unsigned int num_params, unsigned int lw_size, unsigned int rw_size);
  /**
   * @brief y_proj_pla is the interval projection approximation
   * @param s is the series to compress
   * @param num_params is the target dimension for the compression
   * @param lw_size is the size of the first window
   * @param rw_size is the size of the second window
   * @return the Adaptive PLA representation
   */
  Seqddt y_proj_pla(const Seqd& s, unsigned int num_params, unsigned int lw_size, unsigned int rw_size);
}

#endif
