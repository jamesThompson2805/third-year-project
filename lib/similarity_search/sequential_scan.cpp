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

vector<unsigned int> seq_scan::find_k_closest_indexes(const std::vector<double> &series, const std::vector<double> &query, unsigned int k)
{
  if (series.size() <= query.size()) return {0};

  vector<unsigned int> k_closest;
  while (k_closest.size() < k) {
    double min_dist = 1e30;
    unsigned int min_i = -1;
    for (int i=0; i < series.size() - query.size() + 1; ++i) {
      if ( double dist = l2_sqr(series.data() + i, query.data(), query.size()); dist < min_dist ) {
	bool in_closest = false;
	for (unsigned int& k_i : k_closest){ if (k_i == i) in_closest = true; }
	if (!in_closest) k_closest.push_back(i);
      }
    }
    if (min_i == -1) break;
  }
  return k_closest;
}
