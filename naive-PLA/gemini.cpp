#include "./naive_pla.cpp"

/*
 * Gemini Framework to answer Simularity Search
 *
 * For sequence S, apply sliding window size w, apply PLA to each segment to map to 2D feature space.
 * We reason that as the window slides, adjacent subsequences will have similar positions in feature space
 * Hence we form the Minimum Bounding Rectangles around adjacent subsequences, storing start and end positions.
 * Then apply R*-Tree for fast retrieval of MBR's and hence the relevant areas of the sequence to look adjacent
*/

typedef std::vector<float> Sequence;

FloatPair efficienterRegression(const float* const start, const float* const end)
{
  if (end-start<=0) return {0.0, 0.0};  
  int size = end-start+1;

  float x_mean = (size-1.0) / 2.0; 
  float y_mean = std::accumulate(start, end, 0) / (float) size;

  auto calc_residual = [&y_mean](float yi) { return yi - y_mean; };

  float b=0;
  float sqr_x_res=0;
  for (int i=0; i<size; ++i) {
    b += (i - x_mean) * calc_residual( *(start+i) );
    sqr_x_res += (i-x_mean)*(i-x_mean);
  }
  b /= sqr_x_res;
  float a = y_mean - (b*x_mean);

  return {a,b};
}

vector<FloatPair> mapToFeatureSpace(const Sequence& seq, int w)
{


}



int main()
{


  return 0;
}
