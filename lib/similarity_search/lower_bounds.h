#ifndef LOWER_BOUNDS_H
#define LOWER_BOUNDS_H

#include <vector>
#include <tuple>
#include <array>

#include "pla.h"

/**
 * @file lower_bounds.h
 * @brief Header file for distance functions for the indexing scheme created.
 */

typedef std::vector<std::tuple<DoublePair, unsigned int>> Seqddt;
typedef std::vector<double> Seqd;
struct Region {
  DoublePair min_dp;
  unsigned int min_i;
  DoublePair max_dp;
  unsigned int max_i;
};
typedef std::vector<Region> PlaMbr;

/**
 * @brief lower_bounds handles distances to a MBR (Partition Cover) or individual approximation
 */
namespace lower_bounds {
  /**
   * @brief dist_pla_lb_sqr returns the squared l2 distance between series q and compressed s
   * @param q is uncompressed series
   * @param s is compressed series
   */
  double dist_pla_lb_sqr(const Seqd& q, const Seqddt& s);
  /**
   * @brief dist_pla_lb returns the l2 distance between series q and compressed s
   * @param q is uncompressed series
   * @param s is compressed series
   * This is equivalent to std::sqrt( dist_pla_lb_sqr(q,s) )
   */
  double dist_pla_lb(const Seqd& q, const Seqddt& s);
  /**
   * @brief dist_mbr_lb_sqr returns the squared l2 distance between series q and Partition Cover r
   * @param q is uncompressed series
   * @param r is Partition cover
   */
  double dist_mbr_lb_sqr(const Seqd& q, const PlaMbr& r);
  /**
   * @brief dist_mbr_lb returns the l2 distance between series q and Partition Cover r
   * @param q is uncompressed series
   * @param r is Partition cover
   */
  double dist_mbr_lb(const Seqd& q, const PlaMbr& r);

};

#endif
