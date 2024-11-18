#include "double_window.h"
#include "pla.h"

using std::vector;
using std::tuple;

#include <limits>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <array>

#include <iostream>

using std::array;

inline double mean_diff(const double* const d1, const double* const d2, unsigned int lw_size, unsigned int rw_size)
{
  if (d2-d1+1 != lw_size+rw_size || lw_size == 0 || rw_size == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return std::abs( std::accumulate(d1+lw_size, d2+1, 0.0) - std::accumulate(d1, d1 + lw_size, 0.0) );
}

unsigned int greatest_diff_index(const double* const d1, const double* const d2, unsigned int lw_size, unsigned int rw_size)
{
  if (d2-d1+1 < lw_size+rw_size || lw_size == 0 || rw_size == 0) return 0;
  double max_diff = 0.0;
  unsigned int max_diff_index = 0;
  for (int start=0; start <= d2-d1+1 - lw_size - rw_size; ++start) {
    if ( double curr_diff = mean_diff(d1+start, d1+start+lw_size+rw_size-1, lw_size, rw_size); !std::isnan(curr_diff) && curr_diff > max_diff ) { 
      max_diff = curr_diff;
      max_diff_index = start;
    } else if ( std::isnan(curr_diff)) {
      continue;
    }
  }

  return max_diff_index;
}


vector<tuple<double, unsigned int>> d_w::simple_paa(const vector<double> &s, unsigned int num_params, unsigned int lw_size, unsigned int rw_size)
{

  if (lw_size == 0 || rw_size == 0 || num_params <= 1) return {};
  if (s.size() == 0 ) return {};
  unsigned int num_ints = num_params/2;
  
  vector<array<unsigned int,2>> splits;
  splits.push_back( {0, (unsigned int) s.size() - 1} );

  double max_diff = 0.0;
  unsigned int max_diff_split = 0;
  unsigned int max_diff_index = 0;
  unsigned int curr_diff_index = 0;
  while (splits.size() < num_ints) {
    // find location of greatest difference in mean for each split
    // and find maximum difference overall
    max_diff = 0.0;
    max_diff_split = 0;
    max_diff_index = 0;
    curr_diff_index = 0;
    for (int i=0; i<splits.size(); ++i) {
      if ( splits[i][1] - splits[i][0]+1 < lw_size + rw_size) continue; // skip if the interval is too small to examine with the window
      curr_diff_index = greatest_diff_index(s.data() + splits[i][0], s.data() + splits[i][1], lw_size, rw_size);
      if ( double curr_diff = mean_diff(s.data() + splits[i][0] + curr_diff_index, s.data() + splits[i][0] + curr_diff_index + lw_size + rw_size - 1, lw_size, rw_size);
	    curr_diff > max_diff) {
	max_diff = curr_diff;
	max_diff_split = i;
	max_diff_index = curr_diff_index;
      }
    }

    // create a split at that location
    auto [l, r] = splits[max_diff_split];
    splits.erase( splits.begin() + max_diff_split);
    splits.insert( splits.begin() + max_diff_split, { l + max_diff_index + lw_size, r });
    splits.insert( splits.begin() + max_diff_split, { l, l + max_diff_index + lw_size -1 });
  }
  
  vector<tuple<double, unsigned int>> v_paa;
  for (const auto [l, r] : splits) {
    v_paa.push_back( { std::accumulate(s.data() + l, s.data() + r+1, 0.0) / (double) (r-l+1), r});
   }
  return v_paa;
}


inline double int_expansion(const double* const d1, const double* const d2, unsigned int lw_size, unsigned int rw_size)
{
  if (d2-d1+1 != lw_size+rw_size || lw_size == 0 || rw_size == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const auto [min1, max1] = std::minmax_element(d1, d1+lw_size);
  const auto [min2, max2] = std::minmax_element(d1+lw_size, d2+1);
  return std::max(*max2-*max1, 0.0) + std::max(*min1-*min2, 0.0);
}

unsigned int greatest_int_diff_index(const double* const d1, const double* const d2, unsigned int lw_size, unsigned int rw_size)
{
  if (d2-d1+1 < lw_size+rw_size || lw_size == 0 || rw_size == 0) return 0;
  double max_diff = 0.0;
  unsigned int max_diff_index = 0;
  for (int start=0; start <= d2-d1+1 - lw_size - rw_size; ++start) {
    if ( double curr_diff = int_expansion(d1+start, d1+start+lw_size+rw_size-1, lw_size, rw_size); !std::isnan(curr_diff) && curr_diff > max_diff ) { 
      max_diff = curr_diff;
      max_diff_index = start;
    } else if ( std::isnan(curr_diff)) {
      continue;
    }
  }

  return max_diff_index;
}

vector<tuple<double, unsigned int>> d_w::y_proj_paa(const vector<double> &s, unsigned int num_params, unsigned int lw_size, unsigned int rw_size)
{

  if (lw_size == 0 || rw_size == 0 || num_params <= 1) return {};
  if (s.size() == 0 ) return {};
  unsigned int num_ints = num_params/2;
  
  vector<array<unsigned int,2>> splits;
  splits.push_back( {0, (unsigned int) s.size() - 1} );

  double max_diff = 0.0;
  unsigned int max_diff_split = 0;
  unsigned int max_diff_index = 0;
  unsigned int curr_diff_index = 0;
  while (splits.size() < num_ints) {
    // find location of greatest difference in mean for each split
    // and find maximum difference overall
    max_diff = 0.0;
    max_diff_split = 0;
    max_diff_index = 0;
    curr_diff_index = 0;
    for (int i=0; i<splits.size(); ++i) {
      if ( splits[i][1] - splits[i][0]+1 < lw_size + rw_size) continue; // skip if the interval is too small to examine with the window
      curr_diff_index = greatest_int_diff_index(s.data() + splits[i][0], s.data() + splits[i][1], lw_size, rw_size);
      if ( double curr_diff = int_expansion(s.data() + splits[i][0] + curr_diff_index, s.data() + splits[i][0] + curr_diff_index + lw_size + rw_size - 1, lw_size, rw_size);
	    curr_diff > max_diff) {
	max_diff = curr_diff;
	max_diff_split = i;
	max_diff_index = curr_diff_index;
      }
    }

    // create a split at that location
    auto [l, r] = splits[max_diff_split];
    splits.erase( splits.begin() + max_diff_split);
    splits.insert( splits.begin() + max_diff_split, { l + max_diff_index + lw_size, r });
    splits.insert( splits.begin() + max_diff_split, { l, l + max_diff_index + lw_size -1 });
  }
  
  vector<tuple<double, unsigned int>> v_paa;
  for (const auto [l, r] : splits) {
    v_paa.push_back( { std::accumulate(s.data() + l, s.data() + r+1, 0.0) / (double) (r-l+1), r});
   }
  return v_paa;
}

// for x_0, .., x_n, calculates the difference in mean( dx_1, .., dx_lw_size), mean( dx_lw_size+1, dx_rw_size+lw_size)
inline double mean_differential_diff(const double* const d1, const double* const d2, unsigned int lw_size, unsigned int rw_size)
{
  if (lw_size == 0 || rw_size == 0 || d2-d1 != lw_size+rw_size) {
    return std::numeric_limits<double>::quiet_NaN(); 
  }
  double diff_1 = 0.0;
  double diff_2 = 0.0;
  for (int i=0; i<lw_size; ++i) {
    diff_1 += *(d1+i+1) - *(d1+i);
  }
  for (int i=0; i<rw_size; ++i) {
    diff_2 += *(d1+lw_size+i+1) - *(d1+lw_size+i);
  }
  return std::abs( diff_2 - diff_1);
}

unsigned int greatest_differential_index(const double* const d1, const double* const d2, unsigned int lw_size, unsigned int rw_size)
{
  if (d2-d1 < lw_size+rw_size || lw_size == 0 || rw_size == 0) return 0;
  double max_diff = 0.0;
  unsigned int max_diff_index = 0;
  for (int start=1; start <= d2-d1-rw_size-lw_size+1; ++start) {
    if ( double curr_diff = mean_differential_diff(d1+start-1, d1+start+lw_size+rw_size-1, lw_size, rw_size)
	  ; !std::isnan(curr_diff) && curr_diff > max_diff ) { 
      max_diff = curr_diff;
      max_diff_index = start;
    } else if ( std::isnan(curr_diff)) {
      std::cout << "whoops" << std::endl;
      continue;
    }
  }

  return max_diff_index-1;
}

vector<tuple<DoublePair, unsigned int>> d_w::simple_pla(const vector<double> &s, unsigned int num_params, unsigned int lw_size, unsigned int rw_size)
{

  if (lw_size == 0 || rw_size == 0 || num_params <= 2) return {};
  if (s.size() == 0 ) return {};
  unsigned int num_ints = num_params/3;
  
  vector<array<unsigned int,2>> splits;
  splits.push_back( {0, (unsigned int) s.size() - 1} );

  double max_diff = 0.0;
  unsigned int max_diff_split = 0;
  unsigned int max_diff_index = 0;
  unsigned int curr_diff_index = 0;
  bool all_splits_too_small = false;
  while (splits.size() < num_ints) {
    // find location of greatest difference in mean for each split
    // and find maximum difference overall
    max_diff = 0.0;
    max_diff_split = 0;
    max_diff_index = 0;
    curr_diff_index = 0;
    all_splits_too_small = true;
    for (int i=0; i<splits.size(); ++i) {
      if ( splits[i][1] - splits[i][0] < lw_size + rw_size) continue; // skip if the interval is too small to examine with the window
      all_splits_too_small = false;
      curr_diff_index = greatest_differential_index(s.data() + splits[i][0], s.data() + splits[i][1], lw_size, rw_size);

      std::cout << "	split " << splits[i][0] << " " << splits[i][1] << " best index: " << curr_diff_index << std::endl;
      std::cout << "	its differential diff " << mean_differential_diff(s.data() + splits[i][0] + curr_diff_index, s.data() + splits[i][0] + curr_diff_index + lw_size + rw_size, lw_size, rw_size)
<< std::endl;

      if ( double curr_diff = mean_differential_diff(s.data() + splits[i][0] + curr_diff_index, s.data() + splits[i][0] + curr_diff_index + lw_size + rw_size, lw_size, rw_size)
	  ; !std::isnan(curr_diff) && curr_diff > max_diff) {
	std::cout << "	Best split" << std::endl;
	max_diff = curr_diff;
	max_diff_split = i;
	max_diff_index = curr_diff_index;
      } else if (std::isnan(curr_diff)) {
	std::cout << splits[i][0] << " " << splits[i][1] << std::endl;
      }
    }

    if (all_splits_too_small) break;

    // create a split at that location
    auto [l, r] = splits[max_diff_split];
    std::cout << "current diff index: " << max_diff_index << std::endl;
    std::cout << "splitting: (" << l << "," << r << ") into ";
    std::cout << "(" << l << "," << l+max_diff_index + lw_size - 1 << ") and ";
    std::cout << "(" << l+max_diff_index + lw_size << "," << r << ")" << std::endl;
    splits.erase( splits.begin() + max_diff_split);
    splits.insert( splits.begin() + max_diff_split, { l + max_diff_index + lw_size, r });
    splits.insert( splits.begin() + max_diff_split, { l , l + max_diff_index + lw_size -1 });
  }
  
  vector<tuple<DoublePair, unsigned int>> v_paa;
  for (const auto [l, r] : splits) {
    std::cout << l << " : " << r << std::endl;
    v_paa.push_back( { pla::regression(s.data() + l, s.data() + r), r});
   }
  return v_paa;
}

// for x_0, .., x_n, calculates the difference in interval sizes from proj( dx_1, .., dx_lw_size) added proj( dx_lw_size+1, dx_rw_size+lw_size)
inline double differential_int_difference(const double* const d1, const double* const d2, unsigned int lw_size, unsigned int rw_size)
{
  if (lw_size == 0 || rw_size == 0 || d2-d1 != lw_size+rw_size) {
    return std::numeric_limits<double>::quiet_NaN(); 
  }
  double min1 = 100.0, max1 = -100.0;
  double min2 = 100.0, max2 = -100.0;
  double curr;
  for (int i=0; i<lw_size; ++i) {
    curr = *(d1+i+1) - *(d1+i);
    if (curr > max1) max1=curr;
    if (curr < min1) min1=curr;
    
  }
  for (int i=0; i<rw_size; ++i) {
    curr = *(d1+lw_size+i+1) - *(d1+lw_size+i);
    if (curr > max2) max2=curr;
    if (curr < min2) min2=curr;
  }
  return std::max( max2 - max1, 0.0) + std::max( min1 - min2, 0.0);
}

unsigned int greatest_differential_int_index(const double* const d1, const double* const d2, unsigned int lw_size, unsigned int rw_size)
{
  if (d2-d1 < lw_size+rw_size || lw_size == 0 || rw_size == 0) return 0;
  double max_diff = -1.0;
  unsigned int max_diff_index = 0;
  for (int start=1; start <= d2-d1-rw_size-lw_size+1; ++start) {
    if ( double curr_diff = differential_int_difference(d1+start-1, d1+start+lw_size+rw_size-1, lw_size, rw_size)
	  ; !std::isnan(curr_diff) && curr_diff > max_diff ) { 
      max_diff = curr_diff;
      max_diff_index = start;
    }
  }

  return max_diff_index-1;
}

vector<tuple<DoublePair, unsigned int>> d_w::y_proj_pla(const vector<double> &s, unsigned int num_params, unsigned int lw_size, unsigned int rw_size)
{

  if (lw_size == 0 || rw_size == 0 || num_params <= 2) return {};
  if (s.size() == 0 ) return {};
  unsigned int num_ints = num_params/3;
  
  vector<array<unsigned int,2>> splits;
  splits.push_back( {0, (unsigned int) s.size() - 1} );

  double max_diff = 0.0;
  unsigned int max_diff_split = 0;
  unsigned int max_diff_index = 0;
  unsigned int curr_diff_index = 0;
  bool all_splits_too_small = false;
  while (splits.size() < num_ints) {
    // find location of greatest difference in mean for each split
    // and find maximum difference overall
    max_diff = 0.0;
    max_diff_split = 0;
    max_diff_index = 0;
    curr_diff_index = 0;
    all_splits_too_small = true;
    for (int i=0; i<splits.size(); ++i) {
      if ( splits[i][1] - splits[i][0] < lw_size + rw_size) continue; // skip if the interval is too small to examine with the window
      all_splits_too_small = false;
      curr_diff_index = greatest_differential_int_index(s.data() + splits[i][0], s.data() + splits[i][1], lw_size, rw_size);
      std::cout << splits[i][0] << " " << splits[i][1] << " : " << curr_diff_index << std::endl;

      if ( double curr_diff = differential_int_difference(s.data() + splits[i][0] + curr_diff_index, s.data() + splits[i][0] + curr_diff_index + lw_size + rw_size, lw_size, rw_size)
	  ; !std::isnan(curr_diff) && curr_diff > max_diff) {
	max_diff = curr_diff;
	max_diff_split = i;
	max_diff_index = curr_diff_index;
      } else if (std::isnan(curr_diff)) {
	std::cout << splits[i][0] << " " << splits[i][1] << std::endl;
      }
    }

    if (all_splits_too_small) break;

    // create a split at that location
    auto [l, r] = splits[max_diff_split];
    splits.erase( splits.begin() + max_diff_split);
    splits.insert( splits.begin() + max_diff_split, { l + max_diff_index + lw_size, r });
    splits.insert( splits.begin() + max_diff_split, { l , l + max_diff_index + lw_size -1 });
  }
  
  vector<tuple<DoublePair, unsigned int>> v_paa;
  for (const auto [l, r] : splits) {
    std::cout << l << " : " << r << std::endl;
    v_paa.push_back( { pla::regression(s.data() + l, s.data() + r), r});
   }
  return v_paa;
}

