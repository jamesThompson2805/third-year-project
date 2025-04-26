#ifndef CONV_DOUBLE_WINDOW
#define CONV_DOUBLE_WINDOW

#include "pla.h"

/**
 * @file conv_double_window.h is the optimised version of weighted average value double window
 */

/** 
 * @brief c_d_w namespace holds methods for optimised weighted average value double window
 */
namespace c_d_w {
  /**
   * @brief conv_pla calculates the segment split locations by the difference between the gradients of points weighted by distributions l and r
   * @param s is sequence to compress
   * @param num_params is the target dimension to compress to
   * @param l is an array representing a distribution ie. elements are in range [0,1] and sum to 1, this is the distribution for the first window
   * @param r is a distribution for the second window
   * @return an Adaptive PLA representation
   */
  Seqddt conv_pla(const Seqd& s, unsigned int num_params, const Seqd& l, const Seqd& r);
}

#endif
