#include "../R-Tree/RTreeTemplate/RTree.h"
#include "./naive_pla.cpp"

#include <functional>
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

template <int FEATURE_DIM>
class GeminiNaive {
private:
  RTree<IntPair, float, FEATURE_DIM> m_rtree;
  std::function<float(float, float)> m_dist_func; 

public:
  GeminiNaive(const vector<array<float, FEATURE_DIM>>& series, std::function<float(float, float)> dist_func);

  vector<IntPair> candidate_similar_subseq(const vector<array<float, FEATURE_DIM>>& query, float epsilon);
};

template <int FEATURE_DIM>
GeminiNaive<FEATURE_DIM>::GeminiNaive(const vector<array<float, FEATURE_DIM>>& series, std::function<float(float, float)> dist_func)
  : m_dist_func(dist_func)
{
  for (int i=0; i<series.size(); ++i) {
    const auto& arr = series.at(i);
    m_rtree.Insert(arr, arr, {i,i});
  }
}

template <int FEATURE_DIM>
GeminiNaive<FEATURE_DIM> candidate_similar_subseq(const vector<array<float, FEATURE_DIM>>& query, float epsilon)
{
  
}






/* int main()
{


  return 0;
}
*/
