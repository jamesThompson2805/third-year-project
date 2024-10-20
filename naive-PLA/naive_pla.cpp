#include <iostream>
#include <numeric>
#include <algorithm>

#include <chrono>
using namespace std::chrono;

#include <vector>
#include <string>
using std::vector;

#include <matplot/matplot.h>

#include "./tsv_reader.cpp"
#include "./sequential_scan.cpp"
#include "./random_walk.cpp"

#include <filesystem> // included to read directory names of UCRArchive_2018

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

// w is num items compressed to a linear function
vector<FloatPair> sliding_window_regression( const vector<float> series, unsigned int w)
{
  vector<FloatPair> r_pairs;
  for (int i=0; i<series.size() - w; ++i) {
    r_pairs.emplace_back( regression( series.data() + i, series.data() + i + w ) );
  }
  return r_pairs;
}
void plot_coeffs_of_dataset(std::string dataset_name, const vector<FloatPair>& r_pairs, int max_lines=-1)
{
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
  
  save("img/feature_space/" + dataset_name + "_scatter.jpg");
}

void plot_seq_scan_time(std::string dataset_name, float epsilon)
{
  if (epsilon < 0) return;
  vector<float> series = tsv_reader::parse_tsv("../UCRArchive_2018/"+ dataset_name + "/" + dataset_name + "_TEST.tsv");
  std::cout << series.size() << std::endl;
  vector<double> x;
  vector<double> y;
  for (int i=10; i<=1000; i+=10) {
    vector<float> query(i);
    std::copy(series.begin(), series.begin()+i, query.begin());


    auto begin = high_resolution_clock::now(); // https://www.geeksforgeeks.org/measure-execution-time-function-cpp/
    vector<int> subseqs = seq_scan::find_similar_subseq_indexes(series, query, epsilon);
    auto end = high_resolution_clock::now();
    auto time = duration_cast<milliseconds>( end - begin );

    x.emplace_back( (double) i);
    y.emplace_back( (double) time.count());
  }
  matplot::plot(x,y);
  matplot::title(dataset_name + " plot of sequential scan time results (ms)");
  matplot::xlabel("query size");
  matplot::ylabel("time in ms");
  matplot::save("img/seqscans/" + dataset_name + "_seqscan_timed.png");
}

void plot_random_walk(std::function<float ()> gen_func, std::string name, unsigned int len)
{
  RandomWalk walk(gen_func);
  walk.gen_steps(len);

  const auto& walk_container = walk.get_walk();
  vector<double> x;
  vector<double>y;
  int i=0;
  for (const float& f : walk_container) {
    x.emplace_back(i);
    y.emplace_back((double) f);
    i++;
  }
  matplot::plot(x,y);
  matplot::title(name +" random walk");
  matplot::save("img/walks/" + name + " random_walk.png");
}

int main()
{
  vector<std::string> datasets;
  for (const auto& entry : std::filesystem::directory_iterator("../UCRArchive_2018/")) {
    datasets.emplace_back(entry.path().filename());
  }
  std::sort(datasets.begin(), datasets.end());

  vector<float> series = tsv_reader::parse_tsv("../UCRArchive_2018/"+ datasets[0] + "/" + datasets[0] + "_TEST.tsv", 10000);
  vector<float> query(100);
  std::copy(series.begin(), series.begin() + 100, query.begin());

  UnifFunctor unif;
  NormalFunctor normal;
  CauchyFunctor cauchy;
  plot_random_walk(unif, "Uniform", 1000);
  plot_random_walk(normal, "Normal", 1000);
  plot_random_walk(cauchy, "Cauchy", 1000);
  
  return 0;
}
