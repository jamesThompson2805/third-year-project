#include "conv_double_window.h"

#include "pla.h"

using std::vector;
using std::tuple;

#include <queue>
using std::priority_queue;

#include <algorithm>

inline double score( const vector<double>& s, const vector<double>& l, const vector<double>& r, unsigned int i)
{
  double score = 0;
  for (int j=0; j<l.size(); j++) {
    score -= l[j] * ( s[i-l.size()+1+j] - s[i-l.size()+j] );
  }
  for (int j=0; j<r.size(); j++) {
    score += r[j] * ( s[i+j+1] - s[i+j] );
  }
  return score > 0 ? score : -1 * score;
}

vector<tuple<DoublePair, unsigned int>> c_d_w::conv_pla(const vector<double> &s, unsigned int num_params, const vector<double>& l, const vector<double>& r)
{
  unsigned int ns = num_params / 3; // ns is number of segments
  
  auto cmp = [](tuple<unsigned int, double> l, tuple<unsigned int, double> r) { return std::get<1>(l) > std::get<1>(r); };
  priority_queue<tuple<unsigned int, double>, vector<tuple<unsigned int, double>>, decltype(cmp)> p_q(cmp);
  for (int i=0; i<ns-1; i++) { // fill priority q with low score items
    p_q.push( { 0, -1.0 } );
  }

  // using buffer:
  //  want buffer to act like a simple ring buffer
  //  to append to the end, add 1 to buffer_last (MOD buffer.size())
  //   and overwrite at that location in buffer
  //   make sure to overwrite buffer_last with the location you added to
  //  make sure to update 'middle' pointer by add 1 MOD buffer.size()
  vector<double> buffer( l.size() + r.size() - 1, -1); 
  unsigned int buffer_last = buffer.size() - 1;
  unsigned int middle = buffer_last - std::max( r.size(), l.size()-1 ) - 1;
  double buffer_max_val = -1; 
  unsigned int buffer_max_ind = 0;

  unsigned int start = l.size();
  unsigned int end = s.size() - r.size() - 1;

  for (int i=start; i<=end; i++) {
    buffer_last = (buffer_last+1) % buffer.size(); 
    buffer[buffer_last] = score(s, l, r, i);
    middle = (middle + 1) % buffer.size();
    if (buffer_max_ind == buffer_last) { // best was element we overwrote
      for (unsigned int i=0; i<buffer.size(); i++) {
	if (buffer[i] > buffer_max_val) {
	  buffer_max_val = buffer[i];
	  buffer_max_ind = i;
	}
      }
    } // new value is better than current
    buffer_max_ind = buffer[buffer_last] > buffer_max_val ? buffer_last : buffer_max_ind;
    buffer_max_val = std::max( buffer[buffer_last], buffer_max_val );

    if (buffer[middle] != -1 && buffer_max_ind == middle) { // middle is worst element (so add to pq)
      if (std::get<1>( p_q.top() ) < buffer[middle]) {
	p_q.pop();
	p_q.push( { i - std::max(r.size(), l.size()-1) - 1, buffer[middle] } );
      }
    }

  }
  // at end we ensure that elements that could still be splits are added as well
  while (middle != buffer_last) {
    middle = (middle+1) % buffer.size();
    end++;
    if (buffer[middle] != -1 && buffer_max_ind == middle) {
      if (std::get<1>( p_q.top() ) < buffer[middle]) {
	p_q.pop();
	p_q.push( { end - std::max(r.size(), l.size()-1) - 1, buffer[middle] } );
      }
    }
  }
  vector<unsigned int> split_indexes(ns-1);
  for (int i=0; i<ns-1; i++){
    split_indexes[i] = std::get<0>(p_q.top());
    p_q.pop();
  }
  split_indexes.push_back(s.size() - 1);
  std::sort(split_indexes.begin(), split_indexes.end());

  vector<tuple<DoublePair, unsigned int>> apla(ns);
  for (int i=0; i<ns; i++) {
    unsigned int end_index = split_indexes[i];
    if (i==0) {
      apla[i] = { pla::regression(s.data(), s.data()+end_index), end_index };
    } else {
      unsigned int start_index = std::get<1>( apla[i-1] ) + 1;
      apla[i] = { pla::regression(s.data() + start_index, s.data() + end_index), end_index } ;
    }
  }
  
  return apla;
}
