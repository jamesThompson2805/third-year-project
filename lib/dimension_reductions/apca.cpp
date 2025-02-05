#include "apca.h"

#include "paa.h"

#include <algorithm>
#include <cmath>
#include <queue>

#include <iostream>

using std::vector;
using std::tuple;
using std::priority_queue;

vector<tuple<double, unsigned int>> apca::apca(const std::vector<double> &s, unsigned int num_params)
{
  using std::ceil, std::log2, std::abs, std::copy, std::for_each;
  unsigned int num_segments = num_params / 2;
  
  vector<vector<double>> pairwise_means;
  vector<vector<double>> pairwise_diffs;
  pairwise_means.push_back({});
  pairwise_diffs.push_back({});

  // pad out with zeroes
  pairwise_means[0].resize( 1 << (int) ceil( log2( (double) s.size()) ) );
  copy(s.cbegin(), s.cend(), pairwise_means[0].begin());

  // calculate all pairwise differences and means
  while ( pairwise_means.back().size() != 1 ) {
    pairwise_diffs.push_back({});
    pairwise_means.push_back({});
    int iback = pairwise_means.size() - 1;
    for (int i=0; i<pairwise_means[iback - 1].size(); i+=2) {
      pairwise_means[iback].push_back( (double)(pairwise_means[iback-1][i] + pairwise_means[iback-1][i+1]) / 2.0 );
      pairwise_diffs[iback-1].push_back( (double)(pairwise_means[iback-1][i] - pairwise_means[iback-1][i+1]) / 2.0 );
    }
  }
  for (auto& v : pairwise_diffs) { 
    for (auto& d : v) {
      std::cout << d << " ";
    }
    std::cout << std::endl;
  }

  // grab only largest num_segments no. of coefficients 
  // note that index_f automatically normalises the indexed value
  auto normal = [&pairwise_diffs](const int& v){ return pow( sqrt(2), pairwise_diffs.size() -1 -v); };
  auto index_f = [&pairwise_diffs, &normal](const int& v, const int& i){ return pairwise_diffs[v][i] / normal(v); };
  auto comp = [](const tuple<double,int,int>& a, const tuple<double,int,int>& b){ return abs(std::get<0>(a)) > abs(std::get<0>(b)); };

  priority_queue<tuple<double,int,int>, vector<tuple<double,int,int>>, decltype(comp) > priQ(comp);
  priQ.push( {pairwise_means.back()[0], 0, 0} ); // technically res0_value is also a DWT coefficient
  pairwise_means.back()[0] = 0.0;

  for (int v=0; v<pairwise_diffs.size(); v++) {
    for (int i=0; i<pairwise_diffs[v].size(); i++) {
      if ( abs(index_f(v,i)) >= abs(std::get<0>( priQ.top() )) || priQ.size() < num_segments ) {
	priQ.push( { index_f(v,i), v, i } );
	if (priQ.size() > num_segments)
	  priQ.pop();
      }
      pairwise_diffs[v][i] = 0.0;
    }
  }
  // set values back only to chosen coefficients
  while (priQ.size() > 0) {
    auto [d, v,i] = priQ.top();
    priQ.pop();
    if (v == 0 && i == 0) { // kept first DWT coefficient
      pairwise_means.back()[0] = d;
    } else {
      pairwise_diffs[v][i] = d * normal(v);
    }
  }
  for (auto& v : pairwise_diffs) { 
    for (auto& d : v) {
      std::cout << d << " ";
    }
    std::cout << std::endl;
  }

  // reconstruct approximation
  for (int v=pairwise_means.size()-1; v>0; v--) {
    for (int i=0; i<pairwise_means[v].size(); i++) {
      pairwise_means[v-1][2*i] = pairwise_means[v][i] + pairwise_diffs[v-1][i];
      pairwise_means[v-1][2*i+1] = pairwise_means[v][i] - pairwise_diffs[v-1][i];
    }
  }
  for (auto& v : pairwise_means) { 
    for (auto& d : v) {
      std::cout << d << " ";
    }
    std::cout << std::endl;
  }

  // find segments and form actual approximation
  vector<unsigned int> segments;
  for (unsigned int i=0; i<s.size()-2; i++) {
    if (pairwise_means[0][i] != pairwise_means[0][i+1])
      segments.push_back(i);
  }
  segments.push_back(s.size()-1);

  vector<tuple<double, unsigned int>> apca;
  apca.push_back( { paa::get_mean( s.data(), s.data()+segments[0] ), segments[0] } );
  for (unsigned int segi=1; segi<segments.size(); segi++) {
    apca.push_back(
      { paa::get_mean( s.data()+segments[segi-1]+1, s.data()+segments[segi] ), segments[segi] }
    );
  }

  // merge segments until we have the correct number of segments
  while (apca.size() > num_segments) {
    double min_dist = 10000000.0;
    double min_ind = 0;
    for (int segi=0; segi<apca.size()-1; segi++) {
      if ( double dist = abs(std::get<0>( apca[segi] ) - std::get<0>( apca[segi] )); dist < min_dist ) {
	min_dist = dist;
	min_ind = segi;
      }
    }
    if (min_ind == 0) {
      apca[min_ind] = { paa::get_mean(s.data(), s.data() + std::get<1>(apca[min_ind+1])), std::get<1>(apca[min_ind+1]) };
      apca.erase(apca.begin()+min_ind+1);
    } else {
      apca[min_ind] = { paa::get_mean(s.data() + std::get<1>(apca[min_ind-1])+1, s.data() + std::get<1>(apca[min_ind+1])), std::get<1>(apca[min_ind+1]) };
      apca.erase(apca.begin()+min_ind+1);
    }
  }

  return apca;
}
