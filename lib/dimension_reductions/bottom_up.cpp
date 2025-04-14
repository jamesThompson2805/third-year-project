#include "bottom_up.h"

#include <algorithm>

Seqddt bottom_up::bottom_up(const Seqd& s, double eps, ERROR_F err)
{
  Seqddt apla;
  std::vector<double> merge_cost;
  for (int i=0; i<s.size()/2; i++) {
    DoublePair reg_pair = pla::regression(s.data()+2*i, s.data()+2*i+1);
    apla.push_back( { reg_pair, 2*i+1 } );
  }
  if (s.size() % 2 != 0)
    apla.push_back( {{s.back(), 0}, (unsigned int) s.size()-1} );

  auto calc_merge_cost = [&](unsigned int i){
    unsigned int start_i = i==0 ? 0 : std::get<1>(apla[i-1])+1;
    DoublePair regr_comb = pla::regression( s.data()+start_i, s.data() + std::get<1>(apla[i+1]) );
    return err(s.data()+start_i, regr_comb, std::get<1>(apla[i+1]) - start_i );
  };

  unsigned int start_index = 0;
  for (int i=0; i<apla.size()-1; i++) {
    merge_cost.push_back(calc_merge_cost(i));
  }
  auto get_min = [](const Seqd& s){ return std::distance( s.begin(), std::min_element(s.begin(), s.end()) ); };

  while (merge_cost.size() > 0 && merge_cost[ get_min(merge_cost) ] < eps )
  {
    unsigned int i = get_min(merge_cost);
    unsigned int start_i = i==0 ? 0 : std::get<1>(apla[i-1])+1;
    DoublePair regr_comb = pla::regression(s.data() + start_i, s.data() + std::get<1>(apla[i+1]));
    apla[i] = { regr_comb, std::get<1>(apla[i+1]) };

    apla.erase( std::next(apla.begin(), i+1) );
    merge_cost.erase( std::next(merge_cost.begin(), i+1) );
    if ( i != 0 )
      merge_cost[i-1] = calc_merge_cost(i-1);
    if (i < apla.size()-1 )
      merge_cost[i] = calc_merge_cost(i);
  }

  return apla;
}

double bottom_up::se(const double *const s, const DoublePair & dp, unsigned int len)
{
  if ( len == 0 ) return 0;
  double se = 0;
  for (int i=0; i<len; i++) {
    se += (s[i] - dp[0] - dp[1]*i) * (s[i] - dp[0] - dp[1]*i);
  }
  return se;
}
double bottom_up::maxdev( const double *const s, const DoublePair& dp, unsigned int len)
{
  if ( len == 0 ) return 0;
  double maxdev = 0;
  for (int i=0; i<len; i++) {
    maxdev = std::max( maxdev, std::abs(s[i] - dp[0] - dp[1]*i) );
  }
  return maxdev;
}

Seqddt bottom_up::bottom_up_early_cutoff(const Seqd& s, double eps, ERROR_F err, unsigned int num_seg)
{
  Seqddt apla;
  std::vector<double> merge_cost;
  for (int i=0; i<s.size()/2; i++) {
    DoublePair reg_pair = pla::regression(s.data()+2*i, s.data()+2*i+1);
    apla.push_back( { reg_pair, 2*i+1 } );
  }
  if (s.size() % 2 != 0)
    apla.push_back( {{s.back(), 0}, (unsigned int) s.size()-1} );

  auto calc_merge_cost = [&](unsigned int i){
    unsigned int start_index = i==0 ? 0 : std::get<1>(apla[i-1])+1;
    DoublePair regr_of_adj = pla::regression( s.data()+start_index, s.data() + std::get<1>(apla[i+1]) );
    return err(s.data()+start_index, regr_of_adj, std::get<1>(apla[i+1]) - start_index );
  };

  unsigned int start_index = 0;
  for (int i=0; i<apla.size()-1; i++) {
    merge_cost.push_back(calc_merge_cost(i));
  }
  auto get_min = [](const Seqd& s){ return std::distance( s.begin(), std::min_element(s.begin(), s.end()) ); };

  while (merge_cost.size() > 0 && merge_cost[ get_min(merge_cost) ] < eps && apla.size() < num_seg )
  {
    unsigned int i = get_min(merge_cost);
    unsigned int start_i = i==0 ? 0 : std::get<1>(apla[i-1])+1;
    DoublePair regr_comb = pla::regression(s.data() + start_i, s.data() + std::get<1>(apla[i+1]));
    apla[i] = { regr_comb, std::get<1>(apla[i+1]) };

    apla.erase( std::next(apla.begin(), i+1) );
    merge_cost.erase( std::next(merge_cost.begin(), i+1) );
    if ( i != 0 )
      merge_cost[i-1] = calc_merge_cost(i-1);
    if (i < apla.size()-1 )
      merge_cost[i] = calc_merge_cost(i);
  }

  return apla;
}
