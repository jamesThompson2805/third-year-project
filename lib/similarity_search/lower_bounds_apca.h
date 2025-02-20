#ifndef LOWER_BOUNDS_APCA_H
#define LOWER_BOUNDS_APCA_H

#include <array>
#include <numeric>
#include <vector>
#include <tuple>

template <unsigned int S>
using ApcaPt =  std::array<std::tuple<double,unsigned int>,S>;
typedef std::vector<double> Seqd;

namespace apca_bounds {

  struct Region { 
    double min_v;
    unsigned int min_i;
    double max_v;
    unsigned int max_i;
  };
  template <unsigned int S>
  using ApcaMBR = std::array<Region, S>;

  inline double region_area( const Region& r) { return (r.max_v - r.min_v) * (r.max_i - r.min_i); }
  inline Region region_merge( const Region& r1, const Region& r2)
  {
    return { std::min(r1.min_v,r2.min_v), std::min(r1.min_i,r2.min_i), std::max(r1.max_v,r2.max_v), std::max(r1.max_i,r2.max_i) };
  }
  inline double dist_to_region_sqr( const double& qi, const Region& r, unsigned int global_offset)
  {
    if (qi < r.min_v) return (qi-r.min_v)*(qi-r.min_v);
    if (qi > r.max_v) return (qi-r.max_v)*(qi-r.max_v);
    return 0.0;
  }
  inline double dist_to_regions_sqr( const double& qi, const Region* const rstart, const Region* const rend, unsigned int global_offset)
  {
    auto min_f = [&](const double& min_d, const Region& r){
      if (double rd = dist_to_region_sqr(qi,r,global_offset); min_d < 0 || rd < min_d)
	return rd;
      return min_d;
    };
    return std::accumulate(rstart, rend+1, -1.0, min_f);
  }

  template <unsigned int S>
  double mbr_area(const ApcaMBR<S>& mbr)
  {
    return std::accumulate( mbr.begin(), mbr.end(), 0.0, [](const double& d, const Region& r){ return d+region_area(r); });
  }
  template <unsigned int S>
  ApcaMBR<S> mbr_merge(const ApcaMBR<S>& mbr1, const ApcaMBR<S>& mbr2)
  {
    ApcaMBR<S> ret;
    for (int i=0; i<S; i++) {
      ret[i] = region_merge(mbr1[i], mbr2[i]);
    }
    return ret;
  }
  template <unsigned int S>
  double dist_to_mbr_sqr( const Seqd& q, const ApcaMBR<S>& mbr ) {
    int active_start_i = 0;
    int active_end_i = 0;
    while ( active_end_i < mbr.size() && mbr[active_end_i].min_i == 0 ){
      active_end_i++;
    }
    active_end_i--;

    double dist = 0.0;
    for (int i=0; i<q.size(); i++) {
      dist += DistToRegionsSqr(q, &mbr[0] + active_start_i, &mbr[0] + active_end_i, i); 

      // adjust window of regions
      while ( active_start_i < mbr.size() && mbr[active_start_i].max_i < i ) {
	active_start_i++;
      }
      while( active_end_i < mbr.size() && mbr[active_end_i].min_i < i+1 ){
	active_end_i++;
      }
      active_end_i--;
    }

    return dist;
  }
};

#endif
