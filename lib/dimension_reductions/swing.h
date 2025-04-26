#ifndef SWING_H
#define SWING_H
#include "pla.h"

/**
 * @file swing.h is header file containing implementation of SWING Filter
 * Online Piece-wise Linear Approximation of Numerical Streams with Precision Guarantees
 */


/**
 * @brief swing namespace contains swing approximation
 */
namespace swing {
  /**
   * @brief swing converts series s to Adaptive PLA by use of maximum error epsilon
   * @param s is the series to compress
   * @param epsilon is the maximum error value
   * @return approximation of s in the Adaptive PLA format
   */
  Seqddt swing(const Seqd& s, double epsilon);
  /**
   * @brief swing_compr converts series s to Adaptive PLA and exploits that segments share endpoints to reduce space
   * @param s is the series to compress
   * @param epsilon is the maximum error value
   * @return approximation of s, only the gradient given and not y intercept as this can be determined by end of previous segment
   */
  std::vector<std::tuple<double, unsigned int>> swing_compr(const Seqd& s, double epsilon);

};
#endif

