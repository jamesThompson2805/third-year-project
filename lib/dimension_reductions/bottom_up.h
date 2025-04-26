#ifndef BOTTOM_UP_H
#define BOTTOM_UP_H

#include "pla.h"

/**
 * @file bottom_up.h is a header file containing methods for bottom up approximation
 */

/**
 * @brief bottom_up namespace containts methods for bottom up
 */
namespace bottom_up {
  /**
   * @brief ERROR_F is a function taking a pointer to an array, straight line and length of array and returning a distance
   */
  using ERROR_F = double (*)(const double *const, const DoublePair&, unsigned int);
  /**
   * @brief se is function returning squared error of approximation, se is an example of ERROR_F
   */
  double se( const double *const, const DoublePair&, unsigned int);
  /**
   * @brief maxdev is function returning maximum deviation of approximation, maxdev is an example of ERROR_F
   */
  double maxdev( const double *const, const DoublePair&, unsigned int);
  /**
   * @brief bottom_up function calculates bottom up on a sequence using some measurement to assess error of segments
   * @param s is series
   * @param num_params is target dimension of approximation
   * @param err is method to assess error of linear approximation on a segment
   * @return sorted array of segments, each storing line and endpoint
   */
  Seqddt bottom_up(const Seqd& s, double num_params, ERROR_F err);
  /**
   * @brief bottom_up_early_cutoff function operates same as bottom_up but ends early if it reaches the target dimension
   * @param s is series
   * @param num_params is target dimension of approximation
   * @param err is method to assess error of linear approximation on a segment
   * @param k is target dimension the approximation ends early at if it reaches it
   * @return sorted array of segments, each storing line and endpoint
   */
  Seqddt bottom_up_early_cutoff(const Seqd& s, double num_params, ERROR_F err, unsigned int k);
};




#endif
