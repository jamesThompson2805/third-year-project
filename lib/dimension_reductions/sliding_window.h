#ifndef SLIDING_WINDOW_H
#define SLIDING_WINDOW_H

#include "pla.h"

/**
 * @file sliding_window.h contains the adaptive PLA approximation sliding window (where a window represents a prospective segment)
 */

/**
 * @brief sw namespace holds the sliding window approximation
 */
namespace sw {
  /**
   * @brief sliding_window approximates the series q by segments such that every point of a segment is within epsilon of its linear approximation
   * @param q is the series to compress
   * @param epsilon is the maximum error value of the approximation (maxdev(approximation) <= epsilon)
   * @return the Adaptive PLA representation of the series q
   */
  Seqddt sliding_window(const Seqd& q, double epsilon);
};

#endif
