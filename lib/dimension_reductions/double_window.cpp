#include "double_window.h"
#include "pla.h"

using std::vector;
using std::tuple;

#include <limits>
#include <algorithm>
#include <cmath>
#include <array>


using std::array;

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
  unsigned int max_diff_index = 1;
  for (int start=1; start <= d2-d1-rw_size-lw_size+1; ++start) {
    if ( double curr_diff = mean_differential_diff(d1+start-1, d1+start+lw_size+rw_size-1, lw_size, rw_size)
	  ; !std::isnan(curr_diff) && curr_diff >= max_diff ) { 
      max_diff = curr_diff;
      max_diff_index = start;
    } else if ( std::isnan(curr_diff)) {
      continue;
    }
  }

  return max_diff_index-1;
}

vector<tuple<DoublePair, unsigned int>> d_w::simple_pla(const vector<double> &s, unsigned int num_params, unsigned int lw_size, unsigned int rw_size)
{

  if (lw_size == 0 || rw_size == 0 || num_params <= 2) return {};
  if (s.size() == 0 ) return {};
  num_params = std::min( (unsigned int) s.size(), num_params);
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

      if ( double curr_diff = mean_differential_diff(s.data() + splits[i][0] + curr_diff_index, s.data() + splits[i][0] + curr_diff_index + lw_size + rw_size, lw_size, rw_size)
	  ; !std::isnan(curr_diff) && curr_diff > max_diff) {
	max_diff = curr_diff;
	max_diff_split = i;
	max_diff_index = curr_diff_index;
      } else if (std::isnan(curr_diff)) {
      }
    }

    if (all_splits_too_small) break;

    // create a split at that location
    auto [l, r] = splits[max_diff_split];
    splits.erase( splits.begin() + max_diff_split);
    splits.insert( splits.begin() + max_diff_split, { l + max_diff_index + lw_size, r });
    splits.insert( splits.begin() + max_diff_split, { l , l + max_diff_index + lw_size -1 });
  }
  
  vector<tuple<DoublePair, unsigned int>> v_pla;
  for (const auto [l, r] : splits) {
    v_pla.push_back( { pla::regression(s.data() + l, s.data() + r), r});
   }
  return v_pla;
}

// for x_0, .., x_n, calculates the difference in interval sizes from proj( dx_1, .., dx_lw_size) added proj( dx_lw_size+1, dx_rw_size+lw_size)
inline double differential_int_difference(const double* const d1, const double* const d2, unsigned int lw_size, unsigned int rw_size)
{
  if (lw_size == 0 || rw_size == 0 || d2-d1 != lw_size+rw_size) {
    return std::numeric_limits<double>::quiet_NaN(); 
  }
  double min1 = 1e30, max1 = -1e30;
  double min2 = 1e30, max2 = -1e30;
  double curr;
  for (int i=0; i<lw_size; ++i) {
    curr = d1[i+1] - d1[i];
    if (curr > max1) max1=curr;
    if (curr < min1) min1=curr;
    
  }
  for (int i=0; i<rw_size; ++i) {
    curr = d1[lw_size+i+1] - d1[lw_size+i];
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
  num_params = std::min( (unsigned int) s.size(), num_params);
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
    max_diff = -1e30;
    max_diff_split = 0;
    max_diff_index = 0;
    curr_diff_index = 0;
    all_splits_too_small = true;
    for (int i=0; i<splits.size(); ++i) {
      if ( splits[i][1] - splits[i][0] < lw_size + rw_size) continue; // skip if the interval is too small to examine with the window
      all_splits_too_small = false;
      curr_diff_index = greatest_differential_int_index(s.data() + splits[i][0], s.data() + splits[i][1], lw_size, rw_size);

      if ( double curr_diff = differential_int_difference(s.data() + splits[i][0] + curr_diff_index, s.data() + splits[i][0] + curr_diff_index + lw_size + rw_size, lw_size, rw_size)
	  ; !std::isnan(curr_diff) && curr_diff > max_diff) {
	max_diff = curr_diff;
	max_diff_split = i;
	max_diff_index = curr_diff_index;
      } else if (std::isnan(curr_diff)) {
      }
    }

    if (all_splits_too_small) {
      break;
    }

    // create a split at that location
    auto [l, r] = splits[max_diff_split];
    splits.erase( splits.begin() + max_diff_split);
    splits.insert( splits.begin() + max_diff_split, { l + max_diff_index + lw_size, r });
    splits.insert( splits.begin() + max_diff_split, { l , l + max_diff_index + lw_size -1 });
  }
  
  vector<tuple<DoublePair, unsigned int>> v_pla;
  for (const auto [l, r] : splits) {
    v_pla.push_back( { pla::regression(s.data() + l, s.data() + r), r});
  }
  return v_pla;
}

