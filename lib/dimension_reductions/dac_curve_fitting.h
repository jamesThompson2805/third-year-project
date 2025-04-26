#ifndef DAC_CURVE_FITTING_H
#define DAC_CURVE_FITTING_H


#include "pla.h"

/**
 * @file dac_curve_fitting.h is header file for Divide and Conquer AKA. Top Down or RDP
 */

/**
 * @brief dac_curve_fitting namespace holds function for top down approximation
 */
namespace dac_curve_fitting {
  /**
   * @brief dac_linear is the RDP form of the top down algorithm
   * @param s is the series to compress
   * @param epsilon is the maximum error value
   * @return APLA approximation of s
   */
  Seqddt dac_linear( const Seqd& s, double epsilon);
  /**
   * @brief dac_linear_early_cutoff is the RDP form of the top down algorithm that cuts off early if it reaches the target dimension
   * @param s is the series to compress
   * @param epsilon is the maximum error value
   * @param num_seg is number of segments to cut off at
   * @return APLA approximation of s
   */
  Seqddt dac_linear_early_cutoff( const Seqd& s, double epsilon, unsigned int num_seg);
}

#endif
