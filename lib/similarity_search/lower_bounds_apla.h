#ifndef LOWER_BOUNDS_APLA_H
#define LOWER_BOUNDS_APLA_H

#include <array>
#include <tuple>
#include <numeric>
#include <vector>
#include <cmath>

#include "pla.h"

/**
 * @file lower_bounds_apla.h
 * @brief Header file for creation and merging of Partition Covers
 */

template <unsigned int S>
using AplaPt =  std::array<std::tuple<DoublePair, unsigned int>,S>;
typedef std::vector<double> Seqd;

/**
 * @brief apla_bounds namespace holds the functions for creation and merging of Partition Covers
 */
namespace apla_bounds {

  /**
   * @brief Region is a struct storing information about a segment of the Partition Cover
   */
  struct Region { 
    DoublePair min_dp;
    unsigned int min_i;
    DoublePair max_dp;
    unsigned int max_i;
  };
  template <unsigned int S>
  using AplaMBR = std::array<Region, S>;

  using std::max, std::min;
  /**
   * @brief region_area calculates the sum of the area of the region in the plane
   * @param r is the region to determine this area for
   * @return non-negative number representing area of region in plane
   */
  inline double region_area( const Region& r) {
    double width = r.max_i - r.min_i;
    double min_a = r.min_dp[0];
    double min_b=  r.min_dp[0] + width*r.min_dp[1];
    double max_a = r.max_dp[0];
    double max_b=  r.max_dp[0] + width*r.max_dp[1];
    return (min(max_a,max_b) - max(min_a,min_b)) * width
	   + (max(max_a,max_b) - min(max_a,max_b)) * width * 0.5
	   + (max(min_a,min_b) - min(min_a,min_b)) * width * 0.5;
  }
  /**
   * @brief region_merge calculates the merge of two segments (regions)
   * @param r1 is the first region
   * @param r2 is the second region
   * @return new region that contains the previous two
   */
  inline Region region_merge( const Region& r1, const Region& r2)
  {
    int min_i = min(r1.min_i, r2.min_i);
    int max_i = max(r1.max_i, r2.max_i);
    if (max_i == min_i) return { { min(r1.min_dp[0],r2.min_dp[0]),0}
			       , (unsigned int) min_i
			       , { max(r1.max_dp[0],r2.max_dp[0]),0}
			       , (unsigned int) max_i };

    int min_i_rel_r1 = min_i - r1.min_i;
    int max_i_rel_r1 = max_i - r1.min_i;
    int min_i_rel_r2 = min_i - r2.min_i;
    int max_i_rel_r2 = max_i - r2.min_i;

    double r1_lmax = r1.max_dp[1]*min_i_rel_r1 + r1.max_dp[0];
    double r1_rmax = r1.max_dp[1]*max_i_rel_r1 + r1.max_dp[0];
    double r2_lmax = r2.max_dp[1]*min_i_rel_r2 + r2.max_dp[0];
    double r2_rmax = r2.max_dp[1]*max_i_rel_r2 + r2.max_dp[0];

    double r1_lmin = r1.min_dp[1]*min_i_rel_r1 + r1.min_dp[0];
    double r1_rmin = r1.min_dp[1]*max_i_rel_r1 + r1.min_dp[0];
    double r2_lmin = r2.min_dp[1]*min_i_rel_r2 + r2.min_dp[0];
    double r2_rmin = r2.min_dp[1]*max_i_rel_r2 + r2.min_dp[0];

    double r_lmax = max( max(r1_lmax, r2_lmax), max(r1_lmin, r2_lmin) );
    double r_lmin = min( min(r1_lmax, r2_lmax), min(r1_lmin, r2_lmin) );
    double r_rmax = max( max(r1_rmax, r2_rmax), max(r1_rmin, r2_rmin) );
    double r_rmin = min( min(r1_rmax, r2_rmax), min(r1_rmin, r2_rmin) );
    return { {r_lmin, (r_rmin-r_lmin)/double(max_i-min_i)}, (unsigned int) min_i , {r_lmax, (r_rmax-r_lmax)/double(max_i-min_i)}, (unsigned int) max_i};
  }

  /**
   * @brief dist_to_region_sqr takes a current point (global_offset, qi) and calcuates the distance to the region r
   * @param qi is the value of the uncompressed time series at time global_offset
   * @param r is the region to measure the distance to
   * @param global_offset is the time index of qi from the time series
   * @return non negative number representing distance to the region
   */
  inline double dist_to_region_sqr( const double& qi, const Region& r, unsigned int global_offset)
  {
    double r_min_est = r.min_dp[1] * (global_offset - r.min_i) + r.min_dp[0];
    double r_max_est = r.max_dp[1] * (global_offset - r.min_i) + r.max_dp[0];

    if (qi < r_min_est) return (qi-r_min_est)*(qi-r_min_est);
    if (qi > r_max_est) return (qi-r_max_est)*(qi-r_max_est);
    return 0.0;
  }
  /**
   * @brief dist_to_regions_sqr takes a current point (global_offset, qi) and calcuates the distance to the array of regions 
   * @param qi is the value of the uncompressed time series at time global_offset
   * @param rstart is a constant pointer to an array of regions (a partition cover)
   * @param rend is a constant pointer to an array of regions (a partition cover)
   * @param global_offset is the time index of qi from the time series
   * @return non negative number representing distance to the region
   */
  inline double dist_to_regions_sqr( const double& qi, const Region* const rstart, const Region* const rend, unsigned int global_offset)
  {
    auto min_f = [&](const double& min_d, const Region& r){
      if (double rd = dist_to_region_sqr(qi,r,global_offset); min_d < 0 || rd < min_d)
	return rd;
      return min_d;
    };
    return std::accumulate(rstart, rend+1, -1.0, min_f);
  }

  /**
   * @brief mbr_area calculates the total area of the partition cover in the plane as the sum of area of regions
   * @param mbr is a Partition Cover
   */
  template <unsigned int S>
  double mbr_area(const AplaMBR<S>& mbr)
  {
    return std::accumulate( mbr.begin(), mbr.end(), 0.0, [](const double& d, const Region& r){ return d+region_area(r); });
  }
  /**
   * @brief mbr_merge merges two partition covers by region wise merges
   * @param mbr1 is a Partition Cover
   * @param mbr2 is a Partition Cover
   * @return a Partition Cover that contains both Partition Covers
   */
  template <unsigned int S>
  AplaMBR<S> mbr_merge(const AplaMBR<S>& mbr1, const AplaMBR<S>& mbr2)
  {
    AplaMBR<S> ret;
    for (int i=0; i<S; i++) {
      ret[i] = region_merge(mbr1[i], mbr2[i]);
    }
    return ret;
  }
  /**
   * @brief dist_to_mbr_sqr takes a time series q and Partition Cover mbr and returns the distance between the two
   * @param q is a time series (uncompressed)
   * @param mbr is a constant reference to a Partition Cover
   * @return non negative number representing distance to the mbr
   */
  template <unsigned int S>
  double dist_to_mbr_sqr( const Seqd& q, const AplaMBR<S>& mbr ) {
    int active_start_i = 0;
    int active_end_i = 0;
    while ( active_end_i < mbr.size() && mbr[active_end_i].min_i == 0 ){
      active_end_i++;
    }
    active_end_i--;

    double dist = 0.0;
    for (int i=0; i<q.size(); i++) {
      // adjust window of regions
      while ( active_start_i < mbr.size()-1 && mbr[active_start_i].max_i < i ) {
	active_start_i++;
      }
      while( active_end_i < mbr.size()-1 && mbr[active_end_i+1].min_i <= i ){
	active_end_i++;
      }

      dist += dist_to_regions_sqr(q[i], &mbr[0] + active_start_i, &mbr[0] + active_end_i, i); 

    }

    return dist;
  }
  
  /**
   * @brief ptrs_to_region converts an array of doubles into a region that bounds them
   * @param start points to the array beginning
   * @param end points to the arrays end (final element)
   * @param g_start_i is the starting index of this in the larger time series being approximated
   * @return region that bounds the contiguous array of points
   */
  inline Region ptrs_to_region(const double* const start, const double* const end, unsigned int g_start_i)
  {
    if (start == end) return { {start[0], 0}, g_start_i, {start[0], 0}, g_start_i };
    DoublePair dp = pla::regression(start, end);
    if (start+1 == end) return { {dp[0], dp[1]}, g_start_i, {dp[0], dp[1]}, g_start_i+1 };
    //std::cout << " dp " << dp[0] << " : " << dp[1] << " width : " << end - start <<  std::endl;
    int max_i = 0,  min_i = 0;
    double max_v = -1, min_v = -1;
    for (int i=0; i<=end-start; i++) {
      double dist_to_line_sqr = (dp[1]*i - start[i] +dp[0])*(dp[1]*i - start[i] +dp[0]) / std::abs( dp[1]*dp[1] + 1 ) ;
      //std::cout << "		index : " << g_start_i + i <<" val " << start[i] << " dist to line " << dist_to_line_sqr << std::endl;	
      if (dist_to_line_sqr > max_v && start[i] >= dp[0]+dp[1]*i) {
	max_i = i;
	max_v = dist_to_line_sqr;
      }
      if (dist_to_line_sqr > min_v && start[i] <= dp[0]+dp[1]*i) {
	min_i = i;
	min_v = dist_to_line_sqr;
      }
    }
    //std::cout << "	chose min i : " << min_i << ", max i : " << max_i << std::endl;
    int width = end - start;
    Region ret = { {start[min_i]-dp[1]*min_i, dp[1]}
		 , g_start_i
		 , {start[max_i]-dp[1]*max_i, dp[1]}
		 , g_start_i+width };
    for (int i=0; i<=end-start; i++) {
      //std::cout <<"		index : " << i+g_start_i << " dist to region " << dist_to_region_sqr(start[i], ret, g_start_i+i) << std::endl;
    }
    return ret;
  }
  /**
   * @brief vec_to_mbr takes a series q and a Adaptive PLA algorithm and returns a Partition Cover that covers q using the algorithm
   * @param q is the uncompressed time series to cover
   * @param f is the DRT function to q to an approximation
   * @return Partition Cover that covers q
   */
  template <unsigned int S>
  AplaMBR<S> vec_to_mbr(const std::vector<double>& q, pla::APLA_DRT f)
  {
    std::vector<std::tuple<DoublePair, unsigned int>> apla = f(q,3*S);
    AplaMBR<S> mbr;
    int start_i = 0;
    for ( int apla_i = 0; apla_i < S; apla_i++ ) {
      auto& [dp,end_i] = apla[apla_i];
      mbr[apla_i] = ptrs_to_region(q.data()+start_i, q.data()+end_i, start_i);
      //std::cout << mbr[apla_i].min_dp[0] << " " << mbr[apla_i].min_dp[1] << " " << mbr[apla_i].min_i << " " << mbr[apla_i].max_dp[0]<< " " << mbr[apla_i].max_dp[1]<< " " << mbr[apla_i].max_i << std::endl;
      start_i = end_i+1;
    }
    return mbr;
  }
  /**
   * @brief vec_to_subseq_mbrs takes a series q and a Adaptive PLA algorithm, returning an array of PC that cover the subsequences of q
   * @param q is the uncompressed time series to cover
   * @param subseq_size is the desired size of the subsequence
   * @param f is the DRT function to q to an approximation
   * @return array of partition covers, the ith PC covers the ith subsequence
   */
  template <unsigned int S>
  std::vector<AplaMBR<S>> vec_to_subseq_mbrs( const std::vector<double>& q, unsigned int subseq_size, pla::APLA_DRT f)
  {
    std::vector<AplaMBR<S>> subseqs_compr;
    for (int i=0; i<q.size() - subseq_size; i++) {
      std::vector<double> q_i( q.cbegin()+i, q.cbegin()+i+subseq_size); 
      subseqs_compr.push_back( vec_to_mbr<S>(q_i, f) );
    }
    return subseqs_compr;
  }
};

#endif
