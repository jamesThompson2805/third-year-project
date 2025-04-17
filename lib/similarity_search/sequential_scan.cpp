#include "sequential_scan.h"

using std::vector;


double seq_scan::l2_sqr(const double* const s1start, const double* const s2start, unsigned int len)
{
  double mse = 0;
  for (int i=0; i < len; ++i) {
    mse += ( *(s1start+i) - *(s2start+i) ) * ( *(s1start+i) - *(s2start+i) );
  }
  return mse;
}

vector<unsigned int> seq_scan::find_similar_subseq_indexes(const vector<double>& series, const vector<double>& query, double epsilon)
{
  if (series.size() <= query.size()) return {0};
  if (epsilon < 0) return {};

  vector<unsigned int> similar_subseqs;
  for (int i=0; i < series.size() - query.size() + 1; ++i) {
    if ( epsilon * epsilon >= l2_sqr(series.data() + i, query.data(), query.size()) ) {
      similar_subseqs.emplace_back(i);
    }
  }
  return similar_subseqs;
}

#include <queue>
#include <tuple>
using std::priority_queue, std::tuple;

vector<unsigned int> seq_scan::find_k_closest_indexes(const std::vector<double> &series, const std::vector<double> &query, unsigned int k)
{
  if (series.size() <= query.size()) return {0};

  auto cmp = [](const tuple<unsigned int, double>& a, const tuple<unsigned int, double> b){
    return std::get<1>(a) > std::get<1>(b);
  };
  priority_queue<tuple<unsigned int, double>, vector<tuple<unsigned int, double>>, decltype(cmp)> pri_q(cmp);

  for (int i=0; i < series.size() - query.size() + 1; ++i) {
    double dist = l2_sqr(series.data() + i, query.data(), query.size());
    pri_q.push( { i, dist } );
  }

  vector<unsigned int> k_closest;
  for (int i=0; i<k; i++) {
    auto [seq_i,d] = pri_q.top();
    pri_q.pop();
    k_closest.push_back(seq_i);
  }
  return k_closest;
}
