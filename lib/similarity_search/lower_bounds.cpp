#include "lower_bounds.h"

#include <cmath>
#include <array>

inline double sqr(double a) { return a * a; }

inline double dist_lines_sqr(const DoublePair& dp1, const DoublePair& dp2, unsigned int len)
{
  return len * sqr( dp1[0]-dp2[0] ) + (dp1[0]-dp1[1])*(dp1[1]-dp2[1])*len*(len-1) + sqr( dp1[1]-dp2[1] )*len*(len-1)*(2*len-1)/6;
}

double lower_bounds::dist_pla_lb_sqr(const Seqd &q, const Seqddt &s)
{
  double error;
  int start_i = 0;
  for (const auto& [dp,end_i] : s) {
    DoublePair q_reg = pla::regression( q.data()+start_i, q.data()+end_i );
    error += dist_lines_sqr(dp, q_reg, end_i-start_i+1);
    start_i = end_i+1;
  }
  return error;
}

double lower_bounds::dist_pla_lb(const Seqd &q, const Seqddt &s)
{
  return std::sqrt( lower_bounds::dist_pla_lb_sqr(q,s) );
}

double dist_region_lb_sqr_i(const Seqd& q, const Region& r, unsigned int qi)
{
  if (r.min_dp[1] == r.max_dp[1]) { // case gradients same => only contains a line
    if (r.min_dp[1] < 0) {
      return sqr( r.max_dp[0] + (qi-r.min_i)*r.min_dp[1] - q[qi] );
    } else {
      return sqr( r.min_dp[0] + (qi-r.min_i)*r.min_dp[1] - q[qi] );
    }
  } else if (r.min_dp[1] < 0 && r.max_dp[1] > 0) { 
    return std::min( sqr(q[qi]-r.min_dp[0]), sqr(q[qi]-r.max_dp[0]) );
  } else if (r.min_dp[1] >= 0 && r.max_dp[1] > 0) {
    if (double dist1 = q[qi] - r.min_dp[0] + r.min_dp[1]*qi; dist1 < 0) return sqr(dist1);
    if (double dist2 = r.max_dp[0] + r.min_dp[1]*(qi-r.max_i+r.min_i) - q[qi] ; dist2 < 0) return sqr(dist2);
    return 0.0;
  } else if (r.min_dp[1] < 0 && r.max_dp[1] >= 0) {
    if (double dist1 = q[qi] - r.min_dp[0] + r.min_dp[1]*(qi-r.max_i+r.min_i); dist1 < 0) return sqr(dist1);
    if (double dist2 = r.max_dp[0] + r.min_dp[1]*qi - q[qi] ; dist2 < 0) return sqr(dist2);
    return 0.0;
  } else {
    return 0.0;
  }
}

double dist_regions_lb_sqr_i(const Seqd& q, const Region *const r_start, const Region *const r_end, unsigned int qi)
{
  if (r_start > r_end) return -1;
  double min_dist = -1.0;
  for (int i=0; i<=r_end-r_start; i++) {
    double dist = dist_region_lb_sqr_i(q, *(r_start+i), qi);
    if ( min_dist < 0 || dist < min_dist) {
      min_dist = dist;
    }
  }
  return min_dist;
}

double lower_bounds::dist_mbr_lb_sqr(const Seqd& q, const PlaMbr& mbr)
{
  int active_start_i = 0;
  int active_end_i = 0;
  while ( active_end_i < mbr.size() && mbr[active_end_i].min_i == 0 ){
    active_end_i++;
  }
  active_end_i--;

  double dist = 0.0;
  for (int i=0; i<q.size(); i++) {
    dist += dist_regions_lb_sqr_i(q, mbr.data() + active_start_i, mbr.data() + active_end_i, i); 

    // adjust window of regions
    while ( active_start_i < mbr.size() && mbr[active_start_i].max_i < i+1 ) {
      active_start_i++;
    }
    while( active_end_i < mbr.size() && mbr[active_end_i].min_i < i+1 ){
      active_end_i++;
    }
    active_end_i--;
  }

  return dist;
}
