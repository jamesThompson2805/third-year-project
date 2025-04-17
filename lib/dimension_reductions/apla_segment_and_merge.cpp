#include "apla_segment_and_merge.h"

using std::vector;
#include <cmath>

double l2_sqr_error(const double *const arr, const DoublePair &regr, unsigned int length)
{
  double sum = 0;
  for (int i=0; i<length; i++) {
    sum += (arr[i] - (regr[0] + regr[1]*i)) * (arr[i] - (regr[0] + regr[1]*i));
  }
  return sum;
}
double maxdev_error(const double *const arr, const DoublePair &regr, unsigned int length)
{
  double max_err = 0;
  for (int i=0; i<length; i++) {
    max_err = std::max( max_err, std::abs( arr[i] - (regr[0] + regr[1]*i) ) );
  }
  return max_err;
}

void segmerge::merge_1(const Seqd &s, Seqddt &s_compr)
{
  double min_error_incr = 1e20;
  unsigned int min_error_index = 0;
  for (int i=0; i<s_compr.size()-1; i++) {
    unsigned int l_start = i==0 ? 0 : std::get<1>(s_compr[i-1])+1;    
    unsigned int l_end = std::get<1>(s_compr[i]);
    unsigned int r_start = l_end+1;
    unsigned int r_end = std::get<1>(s_compr[i+1]);

    double pair_error = l2_sqr_error(s.data()+l_start, std::get<0>(s_compr[i]), l_end-l_start+1)
			+ l2_sqr_error(s.data()+r_start, std::get<0>(s_compr[i+1]), r_end-r_start+1);
    DoublePair combined = pla::regression(s.data()+l_start, s.data()+r_end);
    double comb_error = l2_sqr_error(s.data()+l_start, combined, r_end-l_start+1);

    if (comb_error - pair_error < min_error_incr) {
      min_error_incr = comb_error - pair_error;
      min_error_index = i;
    }
  }
  unsigned int i = min_error_index;
  unsigned int l_start = i==0 ? 0 : std::get<1>(s_compr[i-1])+1;    
  unsigned int l_end = std::get<1>(s_compr[i]);
  unsigned int r_start = l_end+1;
  unsigned int r_end = std::get<1>(s_compr[i+1]);
  DoublePair combined = pla::regression(s.data()+l_start, s.data()+r_end);

  s_compr.erase( std::next(s_compr.begin(), i) );
  s_compr[i] = { combined, r_end };

}

void segmerge::merge_k(const Seqd &s, Seqddt &s_compr, unsigned int k)
{
  if (k == 0) return;
  for (int i=0; i<k; i++)
    segmerge::merge_1(s,s_compr);
}
void segmerge::merge_to_dim(const Seqd &s, Seqddt &s_compr, unsigned int k)
{
  while (s_compr.size() > k/3)
    segmerge::merge_1(s,s_compr);

  for (int i=0; i<s_compr.size(); i++) {
    unsigned int start = i==0 ? 0 : std::get<1>(s_compr[i-1])+1;    
    unsigned int end = std::get<1>(s_compr[i]);
    DoublePair line = pla::regression(s.data()+start, s.data()+end);
    s_compr[i] = { line, end };
  }

}

void segmerge::segment_1(const Seqd &s, Seqddt &s_compr)
{
  double max_error_loss = -1e20;
  unsigned int split_i = 0;
  unsigned int split_loc = 0;

  for (int i=0; i<s_compr.size(); i++) {
    unsigned int l_start = i==0 ? 0 : std::get<1>(s_compr[i-1])+1;    
    unsigned int r_end = std::get<1>(s_compr[i]);
    if (r_end <= l_start+2) continue;

    double max_error_loss_i = -1e20;
    unsigned int split_loc_i = i;

    for (int j=l_start+1; j<r_end-1; j++) {
      unsigned int l_end = j;
      unsigned int r_start = j+1;

      DoublePair p1 = pla::regression(s.data()+l_start, s.data()+l_end);
      DoublePair p2 = pla::regression(s.data()+r_start, s.data()+r_end);
      double pair_error = l2_sqr_error(s.data()+l_start, p1, l_end-l_start+1)
			  + l2_sqr_error(s.data()+r_start, p2, r_end-r_start+1);
      double comb_error = l2_sqr_error(s.data()+l_start, std::get<0>(s_compr[i]), r_end-l_start+1);

      if (pair_error - comb_error > max_error_loss_i) {
	max_error_loss_i = pair_error - comb_error;
	split_loc_i = j;
      }
    }

    if (max_error_loss_i > max_error_loss) {
      max_error_loss = max_error_loss_i;
      split_loc = split_loc_i;
      split_i = i;
    }
  }
  unsigned int i = split_i;
  unsigned int l_start = i==0 ? 0 : std::get<1>(s_compr[i-1])+1;    
  unsigned int l_end = split_loc;
  unsigned int r_start = split_loc+1;
  unsigned int r_end = std::get<1>(s_compr[i]);
  DoublePair p1 = pla::regression(s.data()+l_start, s.data()+l_end);
  DoublePair p2 = pla::regression(s.data()+r_start, s.data()+r_end);

  s_compr[i] = { p1, l_end };
  s_compr.insert( std::next(s_compr.begin(), i+1), { p2, r_end });
}

void segmerge::segment_k(const Seqd &s, Seqddt &s_compr, unsigned int k)
{
  if (k == 0) return;
  for (int i=0; i<k; i++)
    segment_1(s, s_compr);
}
void segmerge::segment_to_dim(const Seqd &s, Seqddt &s_compr, unsigned int k)
{
  while (s_compr.size() < k/3)
    segmerge::segment_1(s, s_compr);

  for (int i=0; i<s_compr.size(); i++) {
    unsigned int start = i==0 ? 0 : std::get<1>(s_compr[i-1])+1;    
    unsigned int end = std::get<1>(s_compr[i]);
    DoublePair line = pla::regression(s.data()+start, s.data()+end);
    s_compr[i] = { line, end };
  }
}

#include <queue>
using std::priority_queue, std::tuple;
void segmerge::segment_k_opt(const Seqd &s, Seqddt &s_compr, unsigned int k)
{
  if (k==0 || k >= s.size()) return; 

  auto cmp = [&](tuple<unsigned int, unsigned int, double>& a, tuple<unsigned int, unsigned int, double> b){ return std::get<2>(a) < std::get<2>(b); };
  priority_queue<tuple<unsigned int, unsigned int, double>, vector<tuple<unsigned int, unsigned int, double>>, decltype(cmp)> pq(cmp);

  auto calc_seg_split_error_loss = [&](unsigned int seg_i) {
    unsigned int l_start = seg_i==0 ? 0 : std::get<1>(s_compr[seg_i-1])+1;    
    unsigned int r_end = std::get<1>(s_compr[seg_i]);
    if (r_end <= l_start+2) return tuple<unsigned int, double>({0,-1e20});

    double max_error_loss_i = -1e20;
    unsigned int split_loc_i = seg_i;

    for (int j=l_start+1; j<r_end-1; j++) {
      unsigned int l_end = j;
      unsigned int r_start = j+1;

      DoublePair p1 = pla::regression(s.data()+l_start, s.data()+l_end);
      DoublePair p2 = pla::regression(s.data()+r_start, s.data()+r_end);
      double pair_error = l2_sqr_error(s.data()+l_start, p1, l_end-l_start+1)
			  + l2_sqr_error(s.data()+r_start, p2, r_end-r_start+1);
      double comb_error = l2_sqr_error(s.data()+l_start, std::get<0>(s_compr[seg_i]), r_end-l_start+1);

      if (pair_error - comb_error > max_error_loss_i) {
	max_error_loss_i = pair_error - comb_error;
	split_loc_i = j;
      }
    }
    return tuple<unsigned int, double>({split_loc_i, max_error_loss_i});
  };

  for (int seg_i=0; seg_i<s_compr.size(); seg_i++) {
    const auto& [split_i, error_loss] = calc_seg_split_error_loss(seg_i);
    pq.push({seg_i, split_i, error_loss});
  }

  while (s_compr.size() < k) {
    const auto& [seg_i, split_i, e] = pq.top();
    pq.pop();
    unsigned int l_start = seg_i==0 ? 0 : std::get<1>(s_compr[seg_i-1])+1;    
    unsigned int l_end = split_i;
    unsigned int r_start = split_i+1;
    unsigned int r_end = std::get<1>(s_compr[seg_i]);
    DoublePair p1 = pla::regression(s.data()+l_start, s.data()+l_end);
    DoublePair p2 = pla::regression(s.data()+r_start, s.data()+r_end);

    s_compr[seg_i] = { p1, l_end };
    s_compr.insert( std::next(s_compr.begin(), seg_i+1), { p2, r_end });

    const auto& [first_split_i, first_error_loss] = calc_seg_split_error_loss(seg_i);
    pq.push({seg_i, first_split_i, first_error_loss});
    
    const auto& [second_split_i, second_error_loss] = calc_seg_split_error_loss(seg_i+1);
    pq.push({seg_i+1, second_split_i, second_error_loss});
  }
  
}
