#include "bottom_up.h"

#include <algorithm>

#include <iostream>

#include <list>
#include <iterator>
using std::list;

Seqddt bottom_up::bottom_up(const Seqd& s, double eps, ERROR_F err)
{
  //std::cout << "began" << std::endl;

  list<std::tuple<DoublePair, unsigned int>> apla;
  list<double> merge_costs;

  for (int i=0; i<s.size()/2; i++) {
    DoublePair reg_pair = pla::regression(s.data()+2*i, s.data()+2*i+1);
    apla.push_back( { reg_pair, 2*i+1 } );
  }
  if (s.size() % 2 != 0)
    apla.push_back( {{s.back(), 0}, (unsigned int) s.size()-1} );

  auto calc_merge_cost = [&](list<std::tuple<DoublePair, unsigned int>>::iterator it){
    unsigned int start_i;
    if( it==apla.begin()) {
      start_i = 0;
    } else {
      --it;
      start_i = std::get<1>( *it ) + 1;
      ++it;
    }

    ++it;
    unsigned int end_i = std::get<1>( *it );
    --it;
    DoublePair regr_comb = pla::regression( s.data()+start_i, s.data()+end_i );
    return err(s.data()+start_i, regr_comb, end_i - start_i );
  };

  double min_merge_cost = 1e30;
  auto min_merge_it = apla.begin();

  unsigned int start_index = 0;
  int i=0;

  auto last_it = std::prev(apla.end());
  for (auto it=apla.begin(); it!=last_it; ++it) {
    double merge_cost_i = calc_merge_cost(it);
    merge_costs.push_back( merge_cost_i );

    if (min_merge_cost > merge_cost_i) {
      min_merge_cost = merge_cost_i;
      min_merge_it = it;
    }

    i++;
  }

  //std::cout << "got here" << std::endl;

  auto get_min = [](const list<double>& s){ return std::distance( s.begin(), std::min_element(s.begin(), s.end()) ); };

  while (merge_costs.size() > 0 && *std::min_element(merge_costs.begin(), merge_costs.end()) < eps )
  {

    auto min_merge_cost_it = std::min_element(merge_costs.begin(), merge_costs.end());
    min_merge_it = std::next( apla.begin(), std::distance(merge_costs.begin(), min_merge_cost_it) );

    //std::cout << "finding start and end" << std::endl;

    unsigned int start_i;
    if( min_merge_it==apla.begin()) {
      start_i = 0;
    } else {
      --min_merge_it;
      start_i = std::get<1>( *min_merge_it ) + 1;
      ++min_merge_it;
    }
    ++min_merge_it;
    unsigned int end_i = std::get<1>( *min_merge_it );
    --min_merge_it;

    //std::cout << "found start and end " << start_i << " " << end_i << std::endl;

    DoublePair regr_comb = pla::regression(s.data() + start_i, s.data() + end_i);
    *min_merge_it = { regr_comb, end_i };

    //std::cout << "set aplas segment to merge to have correct vals" << std::endl;

    //std::cout << std::get<1>(*min_merge_it) << std::endl;
    //std::cout << *min_merge_cost_it << std::endl;
    //std::cout << "iterators valid" << std::endl;


    apla.erase( std::next(min_merge_it) );

    bool is_first_cost = min_merge_cost_it == merge_costs.begin();
    bool is_last_cost = min_merge_cost_it == std::prev(merge_costs.end());
    auto merge_cost_it_prev = is_first_cost ? merge_costs.begin() : std::prev(min_merge_cost_it);
    merge_costs.erase( min_merge_cost_it );
    merge_cost_it_prev = is_first_cost ? merge_costs.begin() : merge_cost_it_prev;

    //std::cout << "removed the elements after the segment to merge" << std::endl;

    if (is_first_cost) {
      if (merge_costs.size() == 0) break;
      *merge_cost_it_prev = calc_merge_cost(min_merge_it);
    } else if (is_last_cost) {
      *merge_cost_it_prev = calc_merge_cost( std::prev(min_merge_it) );
    } else {
      *merge_cost_it_prev = calc_merge_cost( std::prev(min_merge_it) );
      *std::next(merge_cost_it_prev) = calc_merge_cost(min_merge_it );
    }

  }

  Seqddt apla_vec( apla.begin(), apla.end() );
  return apla_vec;
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

  while (merge_cost.size() > 0 && merge_cost[ get_min(merge_cost) ] < eps)
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
