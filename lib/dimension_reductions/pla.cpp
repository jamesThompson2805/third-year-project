#include "pla.h"

#include <numeric>

using std::vector;

vector<vector<float>> pla::chunk_series(vector<float> series, unsigned int chunk_size)
{
  vector<vector<float>> chunked;
  vector<float> chunk;
  int num_chunks = (series.size() / chunk_size) + (series.size()%chunk_size>0);
  for (int i=0; i<num_chunks; ++i) {
    chunk.clear();
    for (int j=0; j<chunk_size; ++j) {
      if (i*chunk_size + j >= series.size()) break;
      chunk.emplace_back( series[i*chunk_size + j] );
    }
    chunked.emplace_back(chunk);
  }
  return chunked;
}

FloatPair pla::regression(const float* const start, const float* const end)
{
  if (end-start<0) return {0.0, 0.0};  
  if (end-start==0) return {*start, 0.0};
  int size = end-start+1;

  float x_mean = (size-1.0) / 2.0; 
  float y_mean = std::accumulate(start, end+1, 0) / (float) size;

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

// w is num items compressed to a linear function
vector<FloatPair> pla::sliding_window_regression( const vector<float> series, unsigned int w)
{
  vector<FloatPair> r_pairs;
  for (int i=0; i<series.size() - w; ++i) {
    r_pairs.emplace_back( regression( series.data() + i, series.data() + i + w ) );
  }
  return r_pairs;
}
