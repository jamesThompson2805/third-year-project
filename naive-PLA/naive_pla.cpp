#include <iostream>
#include <numeric>
#include <algorithm>

#include <vector>
#include <string>
using std::vector;

#include <matplot/matplot.h>

#include "tsv_reader.cpp"

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

FloatPair regression(const float* const start, const float* const end)
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

void plot_coeffs_of_dataset(std::string dataset_name, int max_lines=-1)
{
  vector<float> series = tsv_reader::parse_tsv("../UCRArchive_2018/"
						+ dataset_name
						+ "/"
						+ dataset_name
						+ "_TEST.tsv", max_lines);

  vector<FloatPair> r_pairs;
  int w = 5; // w is num items compressed to a linear function
  for (int i=0; i<series.size() - w; ++i) {
    r_pairs.emplace_back( regression( series.data() + i, series.data() + i + w ) );
  }


  std::cout << r_pairs[0].a << " " << r_pairs[0].b << std::endl;
  std::cout << r_pairs.size() << std::endl;

  using namespace matplot;
  vector<double> x;
  vector<double> y;
  std::for_each(r_pairs.begin(), r_pairs.end(), [&](auto pair){
      x.emplace_back((double) pair.a);
      y.emplace_back((double) pair.b);});

  auto c = linspace(1,10,x.size());
  auto l = scatter(x,y,6,c);
  l->marker_face(true);

  title(dataset_name + "plot of coefficients on (a,b) space a+tb");
  xlabel("constant coefficient a");
  ylabel("first degree coefficient b");
  
  save("img/" + dataset_name + "_scatter.jpg");
}

int main()
{
  vector<std::string> datasets = {"ACSF1", "Adiac", "AllGestureWiimoteX", "ArrowHead", "Beef", "BeetleFly", "BirdChicken", "BME", "Car", "CBF" };

  std::for_each(datasets.begin(), datasets.end(), [](auto name){
      plot_coeffs_of_dataset(name, 1000);
  } );
  return 0;
}
