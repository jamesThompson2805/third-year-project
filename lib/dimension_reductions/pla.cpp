#include "pla.h"

#include <numeric>

using std::vector;

vector<vector<double>> pla::chunk_series(vector<double> series, unsigned int chunk_size)
{
  vector<vector<double>> chunked;
  vector<double> chunk;
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

// for a + tb
DoublePair pla::regression(const double* const start, const double* const end)
{
  if (end-start<0) return {0.0, 0.0};  
  if (end-start==0) return {*start, 0.0};
  int size = end-start+1;

  double x_mean = (size-1.0) / 2.0; 
  double y_mean = std::accumulate(start, end+1, 0.0) / (double) size;

  auto calc_residual = [&y_mean](double yi) { return yi - y_mean; };

  double b=0;
  double sqr_x_res=0;
  for (int i=0; i<size; ++i) {
    b += (i - x_mean) * calc_residual( *(start+i) );
    sqr_x_res += (i-x_mean)*(i-x_mean);
  }
  b =  b / sqr_x_res;
  double a = y_mean - b*x_mean;

  return {a,b};
}
double pla::regression_thru_point(double first_point, const double* const start, const double* const end)
{
  // using formulas from 
  //  https://stats.libretexts.org/Bookshelves/Computing_and_Modeling/Supplemental_Modules_(Computing_and_Modeling)/Regression_Analysis/Simple_linear_regression/Regression_through_the_origin
  if (end-start<0) return 0.0;  
  if (end-start==0) return start[0]-first_point;
  int size = end-start+1;

  double b_hat=0;
  double sqr_x_sum=0;
  for (int i=0; i<size; ++i) {
    b_hat += (i+1) * (start[i]-first_point); // translate each point as to have first point be origin
    sqr_x_sum += (i+1)*(i+1);
  }
  b_hat =  b_hat / sqr_x_sum;
  // no need to translate back as b_hat is the gradient and therefore invariant of location
  return b_hat;
}

// w is num items compressed to a linear function
vector<DoublePair> pla::sliding_window_regression( const vector<double>& series, unsigned int w)
{
  vector<DoublePair> r_pairs;
  for (int i=0; i<series.size() - w; ++i) {
    r_pairs.emplace_back( regression( series.data() + i, series.data() + i + w ) );
  }
  return r_pairs;
}

vector<DoublePair> pla::chunk_regression( const vector<double>& series, unsigned int num_params)
{
  unsigned int interval_size = (2*series.size() / num_params) + ( (2*series.size()) % num_params != 0);
  vector<DoublePair> r_pairs;
  for (int i=0; i<series.size(); i+=interval_size) {
    if (i+interval_size >= series.size()) {
      r_pairs.emplace_back( regression( series.data()+i, series.data()+series.size()-1 ) );
    } else {
      r_pairs.emplace_back( regression( series.data()+i, series.data()+i+interval_size-1 ) );
    }
  }
  return r_pairs;
}

vector<DoublePair> pla::pla(const std::vector<double> &series, unsigned int num_params)
{
  return chunk_regression(series, num_params);
};

double pla::pla_mse(const vector<double> &series, unsigned int num_params)
{
  vector<DoublePair> pla_series = pla::chunk_regression(series, num_params);
  unsigned int interval_size = (2*series.size() / num_params) + ( (2*series.size()) % num_params != 0);

  double mse = 0.0;
  double pla_estimate = 0.0;
  for (int i=0; i< series.size(); ++i) {
    pla_estimate = pla_series[i / interval_size][0] + pla_series[i / interval_size][1] * (i%interval_size);
    mse += (series[i] - pla_estimate) * (series[i] - pla_estimate);
  }
  return mse / series.size();
}

vector<double> pla::pla_to_seq(const std::vector<DoublePair> pla_seq, unsigned int int_size)
{
  vector<double> seq;
  for (const auto& [a,b] : pla_seq) {
    for (int i=0; i<int_size; ++i)
      seq.emplace_back( a + b* (double)i);
  }
  return seq;
}

vector<double> pla::apla_to_seq(const std::vector<std::tuple<DoublePair, unsigned int>> apla_seq)
{
  vector<double> seq;
  int last_index = -1;
  for (const auto& [dp, int_end] : apla_seq) {
    for (int i=0; i+last_index+1 <= int_end; ++i)
      seq.emplace_back(dp[0] + dp[1] * (double)i);
    last_index = int_end;
  }
  return seq;
}
