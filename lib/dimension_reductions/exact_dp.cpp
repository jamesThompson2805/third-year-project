#include "exact_dp.h"

#include <numeric>

using std::vector;
using std::tuple;


inline double quick_mean(const double *const f1, const double *const f2)
{
  double a = 0;
  for (int i=0; i<=f2-f1; i++) {
    a += *(f1+i);
  }
  return a / double (f2 - f1 + 1);
}
inline double quick_se_with_mean( const double *const f1, const double *const f2, double mean)
{
  double se = 0;
  for (int i=0; i<=f2-f1; i++) {
    se += (*(f1+i) - mean) * (*(f1+i) - mean);
  }
  return se;
}

vector< tuple< double, unsigned int>> exact_dp::min_mse_paa( const vector<double>& s, unsigned int num_params)
{
  int num_segments = num_params/2;
  vector< vector< double>> eps( s.size(), vector<double>(num_segments+1,-1.0) );
  vector< vector< double>> lline( s.size(), vector<double>(num_segments+1) );
  vector< vector< double>> first_i( s.size(), vector<double>(num_segments+1) );

  for (int w=0; w<s.size(); w++) {
    double mean = quick_mean( s.data(), s.data()+w);
    lline[w][1] = mean;
    eps[w][1] = quick_se_with_mean( s.data(), s.data()+w, mean);
    first_i[w][1] = 0;
  }

  for (int t=2; t<=num_segments; t++) {
    for (int w=t-1; w<s.size(); w++) {
      for (int a=t-2; a<w; a++) {
	double mean_of_remaining = quick_mean( s.data()+a+1, s.data()+w );
	double err = quick_se_with_mean(s.data()+a+1, s.data()+w, mean_of_remaining) + eps[a][t-1];

	if ( eps[w][t] == -1.0 || err < eps[w][t] ){
	  lline[w][t] = mean_of_remaining;
	  eps[w][t] = err;
	  first_i[w][t] = a+1;
	}

      }
    }
  }

  vector<tuple<double, unsigned int>> paa(num_segments);
  int w = s.size()-1;
  int segment = num_segments-1;
  while (segment >= 0) {
    paa[segment] = { lline[w][segment+1], w };
    w = first_i[w][segment+1] - 1;
    segment--;
  }

  return paa;
}

DoublePair regress(const double* const start, const double* const end)
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

inline double quick_se_with_regress( const double *const f1, const double *const f2, const DoublePair& regressed)
{
  const auto& [a,b] = regressed;
  double se = 0;
  for (int i=0; i<=f2-f1; i++) {
    se += (*(f1+i) - a - b*i) * (*(f1+i) - a - b*i);
  }
  return se;
}

vector< tuple< DoublePair, unsigned int>> exact_dp::min_mse_pla( const vector<double>& s, unsigned int num_params)
{
  int num_segments = num_params/3;
  vector< vector< double>> eps( s.size(), vector<double>(num_segments+1,-1.0) );
  vector< vector<DoublePair>> lline( s.size(), vector<DoublePair>(num_segments+1) );
  vector< vector< double>> first_i( s.size(), vector<double>(num_segments+1) );

  for (int w=0; w<s.size(); w++) {
    DoublePair line = regress( s.data(), s.data()+w);
    lline[w][1] = line;
    eps[w][1] = quick_se_with_regress( s.data(), s.data()+w, line);
    first_i[w][1] = 0;
  }

  for (int t=2; t<=num_segments; t++) {
    for (int w=t-1; w<s.size(); w++) {
      for (int a=t-2; a<w; a++) {
	DoublePair line_of_remain = regress( s.data()+a+1, s.data()+w );
	double err = quick_se_with_regress(s.data()+a+1, s.data()+w, line_of_remain) + eps[a][t-1];

	if ( eps[w][t] == -1.0 || err < eps[w][t] ){
	  lline[w][t] = line_of_remain;
	  eps[w][t] = err;
	  first_i[w][t] = a+1;
	}

      }
    }
  }

  vector<tuple<DoublePair, unsigned int>> pla(num_segments);
  int w = s.size()-1;
  int segment = num_segments-1;
  while (segment >= 0) {
    pla[segment] = { lline[w][segment+1], w };
    w = first_i[w][segment+1] - 1;
    segment--;
  }

  return pla;
}

inline double quick_maxdev_with_mean( const double *const f1, const double *const f2, double mean)
{
  double maxdev = 0;
  for (int i=0; i<=f2-f1; i++) {
    double dev = (*(f1+i) - mean);
    dev = dev >= 0 ? dev : -1.0 * dev;
    maxdev = std::max( dev, maxdev );
  }
  return maxdev;
}
vector< tuple< double, unsigned int>> exact_dp::min_maxdev_paa( const vector<double>& s, unsigned int num_params)
{
  int num_segments = num_params/2;
  vector< vector< double>> eps( s.size(), vector<double>(num_segments+1,-1.0) );
  vector< vector< double>> lline( s.size(), vector<double>(num_segments+1) );
  vector< vector< double>> first_i( s.size(), vector<double>(num_segments+1) );

  for (int w=0; w<s.size(); w++) {
    double mean = quick_mean( s.data(), s.data()+w);
    lline[w][1] = mean;
    eps[w][1] = quick_se_with_mean( s.data(), s.data()+w, mean);
    first_i[w][1] = 0;
  }

  for (int t=2; t<=num_segments; t++) {
    for (int w=t-1; w<s.size(); w++) {
      for (int a=t-2; a<w; a++) {
	double mean_of_remaining = quick_mean( s.data()+a+1, s.data()+w );
	double err = quick_maxdev_with_mean(s.data()+a+1, s.data()+w, mean_of_remaining) + eps[a][t-1];

	if ( eps[w][t] == -1.0 || err < eps[w][t] ){
	  lline[w][t] = mean_of_remaining;
	  eps[w][t] = err;
	  first_i[w][t] = a+1;
	}

      }
    }
  }

  vector<tuple<double, unsigned int>> paa(num_segments);
  int w = s.size()-1;
  int segment = num_segments-1;
  while (segment >= 0) {
    paa[segment] = { lline[w][segment+1], w };
    w = first_i[w][segment+1] - 1;
    segment--;
  }

  return paa;
}

inline double quick_maxdev_with_regress( const double *const f1, const double *const f2, const DoublePair& regressed)
{
  const auto& [a,b] = regressed;
  double maxdev = 0;
  for (int i=0; i<=f2-f1; i++) {
    double dev = (*(f1+i) - a - b*i);
    dev = dev >= 0 ? dev :  -1.0 * dev;
    maxdev = std::max( maxdev, dev );
  }
  return maxdev;
}
vector< tuple< DoublePair, unsigned int>> exact_dp::min_maxdev_pla( const vector<double>& s, unsigned int num_params)
{
  int num_segments = num_params/3;
  vector< vector< double>> eps( s.size(), vector<double>(num_segments+1,-1.0) );
  vector< vector<DoublePair>> lline( s.size(), vector<DoublePair>(num_segments+1) );
  vector< vector< double>> first_i( s.size(), vector<double>(num_segments+1) );

  for (int w=0; w<s.size(); w++) {
    DoublePair line = regress( s.data(), s.data()+w);
    lline[w][1] = line;
    eps[w][1] = quick_maxdev_with_regress( s.data(), s.data()+w, line);
    first_i[w][1] = 0;
  }

  for (int t=2; t<=num_segments; t++) {
    for (int w=t-1; w<s.size(); w++) {
      for (int a=t-2; a<w; a++) {
	DoublePair line_of_remain = regress( s.data()+a+1, s.data()+w );
	double err = quick_se_with_regress(s.data()+a+1, s.data()+w, line_of_remain) + eps[a][t-1];

	if ( eps[w][t] == -1.0 || err < eps[w][t] ){
	  lline[w][t] = line_of_remain;
	  eps[w][t] = err;
	  first_i[w][t] = a+1;
	}

      }
    }
  }

  vector<tuple<DoublePair, unsigned int>> pla(num_segments);
  int w = s.size()-1;
  int segment = num_segments-1;
  while (segment >= 0) {
    pla[segment] = { lline[w][segment+1], w };
    w = first_i[w][segment+1] - 1;
    segment--;
  }

  return pla;
}
