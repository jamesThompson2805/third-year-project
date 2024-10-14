#include <iostream>
#include <numeric>
#include <algorithm>

#include <vector>
using std::vector;


vector<vector<float>> chunk_series(vector<float> series, unsigned int chunk_size)
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

struct FloatPair { float a; float b; };

// Assume series is pairs (0,x0), (1,x1), ..., (n,xn)
FloatPair regression(const vector<float>& series)
{
  if (series.size()==0) return {0.0,0.0};
  float x_mean = float(series.size()-1) / 2.0;
  float y_mean = std::accumulate(series.begin(), series.end(), 0) / (float) series.size();

  auto calc_residual = [&y_mean](float yi) { return yi - y_mean; };

  float b=0;
  float sqr_x_res=0;
  for (int i=0; i<series.size(); ++i) {
    b += (i - x_mean) * calc_residual(series[i]);
    sqr_x_res += (i-x_mean)*(i-x_mean);
  }
  b /= sqr_x_res;
  float a = y_mean - (b*x_mean);

  return {a,b};
}

int main()
{

  vector<float> series = {1, 2, 1, 4, 5, 6, 7, 8, 9, 10};

  vector<vector<float>> chunked = chunk_series(series, 3);
  std::cout << chunked.size() << std::endl;


  vector<FloatPair> r_pairs;
  std::for_each(chunked.begin(), chunked.end(), [&r_pairs](auto interval) {
    r_pairs.emplace_back( regression(interval) );
  });

  std::cout << r_pairs[0].a << " " << r_pairs[0].b << std::endl;

  return 0;
}
