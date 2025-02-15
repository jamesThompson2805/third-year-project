#ifndef LOWER_BOUNDS_H
#define LOWER_BOUNDS_H

#include <vector>
#include <tuple>
#include <array>

#include "pla.h"

typedef std::vector<std::tuple<DoublePair, unsigned int>> Seqddt;
typedef std::vector<double> Seqd;
struct Region {
  DoublePair min_dp;
  unsigned int min_i;
  DoublePair max_dp;
  unsigned int max_i;
};
typedef std::vector<Region> PlaMbr;

namespace lower_bounds {

  double dist_pla_lb_sqr(const Seqd& q, const Seqddt& s);
  double dist_pla_lb(const Seqd& q, const Seqddt& s);
  double dist_mbr_lb_sqr(const Seqd& q, const PlaMbr& r);
  double dist_mbr_lb(const Seqd& q, const PlaMbr& r);

};

#endif
