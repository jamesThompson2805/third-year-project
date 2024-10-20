
#include "../R-Tree/RTreeTemplate/RTree.h"
#include "./naive_pla.cpp"

#include <functional>
#include <algorithm>
#include <array>
using std::array;

/*
 * Gemini Framework to answer Simularity Search
 *
 * For sequence S, apply sliding window size w, apply PLA to each segment to map to 2D feature space.
 * We reason that as the window slides, adjacent subsequences will have similar positions in feature space
 * Hence we form the Minimum Bounding Rectangles around adjacent subsequences, storing start and end positions.
 * Then apply R*-Tree for fast retrieval of MBR's and hence the relevant areas of the sequence to look adjacent
*/

typedef array<int, 2> IntPair;

// Euclidean distance
template <int FEATURE_DIM>
class GeminiNaive {
private:
  RTree<IntPair, float, FEATURE_DIM> m_rtree;
  vector<float> m_series;
  std::function< vector<array<float, FEATURE_DIM>>( vector<float> )> m_dim_reduct;

public:
  GeminiNaive(const vector<array<float, FEATURE_DIM>>& series);

  vector<IntPair> candidate_similar_subseq(const array<float, FEATURE_DIM>& query, float epsilon);
};

template <int FEATURE_DIM>
GeminiNaive<FEATURE_DIM>::GeminiNaive(const vector<array<float, FEATURE_DIM>>& series)
{
  for (int i=0; i<series.size(); ++i) {
    const auto& arr = series.at(i);
    m_rtree.Insert(arr, arr, {i,i});
  }
}

// will use a feature of euclidean space: 
template <int FEATURE_DIM>
vector<IntPair> GeminiNaive<FEATURE_DIM>::candidate_similar_subseq(const array<float, FEATURE_DIM>& query, float epsilon)
{
  if (epsilon < 0.0) return {};
  array<float, FEATURE_DIM> min;
  array<float, FEATURE_DIM> max;
  std::transform(query.cbegin(), query.cend(), min.begin(),[&](const auto& f){return f-epsilon;});
  std::transform(query.cbegin(), query.cend(), max.begin(),[&](const auto& f){return f+epsilon;});
  vector<IntPair> candidates;

  auto retrieve = [&candidates](const IntPair& p){ candidates.emplace_back(candidates); return true;};
  m_rtree.Search(min, max, retrieve);

}






/* int main()
{


  return 0;
}
*/
