#include "sequential_scan.h"

using std::vector;


float seq_scan::mse_l2(const float* const s1start, const float* const s2start, unsigned int len)
{
  float mse = 0;
  for (int i=0; i < len; ++i) {
    mse += ( *(s1start+i) - *(s2start+i) ) * ( *(s1start+i) - *(s2start+i) );
  }
  return mse;
}

vector<int> seq_scan::find_similar_subseq_indexes(const vector<float>& series, const vector<float>& query, float epsilon)
{
  if (series.size() <= query.size()) return {0};
  if (epsilon < 0) return {};

  vector<int> similar_subseqs;
  for (int i=0; i < series.size() - query.size() + 1; ++i) {
    if ( epsilon * epsilon >= mse_l2(series.data() + i, query.data(), query.size()) ) {
      similar_subseqs.emplace_back(i);
    }
  }
  return similar_subseqs;
}
