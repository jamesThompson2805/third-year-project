#include "dac_curve_fitting.h"

#include "pla.h"

#include <numeric>
#include <algorithm>
#include <cmath>
#include <array>

#include <iostream>

using std::vector;
using std::tuple;
using std::array;

double mean(const double* const f1, const double* const f2)
{
  if (f2 < f1) throw ("Not pointing to valid container");
  int size = f2 - f1 + 1;
  return std::accumulate(f1, f2+1, 0.0) / (double) size;
}

vector<tuple<double, unsigned int>> dac_curve_fitting::dac_constant( const vector<double>& series, double epsilon)
{
  if (epsilon < 0) return {};
  if (series.size() == 0) return {};

  vector<array<int, 2>> queue = {};
  queue.push_back( {0, (int) series.size() - 1} );
  vector<tuple<double, unsigned int>> curves;

  while (queue.size() > 0) {
    array<int, 2> intrvl = queue.back();
    queue.pop_back();
    double s_mean = mean( series.data() + intrvl[0], series.data() + intrvl[1] );
    
    const auto max_dev_el = std::max_element(series.begin() + intrvl[0], series.begin() + intrvl[1] + 1,
	[&s_mean](double i, double j){ return std::abs(i-s_mean) < std::abs(j-s_mean); }
    );

    if ( std::abs( *max_dev_el - s_mean) <= epsilon) {
      curves.push_back( {s_mean, intrvl[1] } );
    } else {
      // break_point is relative to the beginning of the subsequence starting series.begin() + intrvl[0]
      int break_point = std::distance(series.begin()+intrvl[0], max_dev_el); 

      if ( break_point == 0) {
	curves.push_back( {*max_dev_el, intrvl[0]} );
	queue.push_back( {intrvl[0]+1,intrvl[1]});
      } else if ( break_point == intrvl[1] - intrvl[0]){
	curves.push_back( {*max_dev_el, intrvl[1]} );
	queue.push_back( {intrvl[0],intrvl[1]-1});
      } else {

	double s1_mean = mean( series.data() + intrvl[0], series.data() + intrvl[0] + break_point - 1);
	double s2_mean = mean( series.data() + intrvl[0] + break_point + 1, series.data() + intrvl[1]);

	bool add_break_to_s1 = std::abs( *max_dev_el - s1_mean) < std::abs(*max_dev_el - s2_mean);
	if (add_break_to_s1) {
	  queue.push_back( {intrvl[0], intrvl[0] + break_point});
	  queue.push_back( {intrvl[0] + break_point + 1,intrvl[1]});
	} else {
	  queue.push_back( {intrvl[0],intrvl[0] + break_point - 1});
	  queue.push_back( {intrvl[0] + break_point,intrvl[1]});
	}
      }
    }
  }
  std::sort(curves.begin(), curves.end(), [](const auto& tp1, const auto& tp2) { return std::get<1>(tp1) < std::get<1>(tp2); });
  return curves;

}

vector<tuple<DoublePair, unsigned int>> dac_curve_fitting::dac_linear( const vector<double>& series, double epsilon)
{
  if (epsilon < 0) return {};
  if (series.size() == 0) return {};

  vector<array<unsigned int, 2>> stack = {};
  stack.push_back( {0, (unsigned int) series.size() - 1} );
  vector<tuple<DoublePair, unsigned int>> curves;

  auto dist_from_line = [](DoublePair line, double x, double y) { return std::abs( y - line[0] - line[1]*x); };

  while (stack.size() > 0) {
    array<unsigned int, 2> intrvl = stack.back();
    stack.pop_back();
    DoublePair regressed = pla::regression( series.data() + intrvl[0], series.data() + intrvl[1] );
    
    unsigned int max_i = intrvl[0];
    double max_dist = 0;
    for (int i=intrvl[0]; i<=intrvl[1]; ++i) {
      if ( dist_from_line(regressed, i-intrvl[0], series[i]) >= dist_from_line( regressed, max_i-intrvl[0], series[max_i])) {
	max_i = i;
	max_dist = dist_from_line(regressed, i-intrvl[0], series[i]);
      }
    }

    if ( max_dist <= epsilon) {
      curves.push_back( {regressed, intrvl[1] } );
    } else {
      if ( max_i == intrvl[0]) {
	curves.push_back( { {series[max_i], 0 }, intrvl[0]} );
	stack.push_back( {intrvl[0]+1,intrvl[1]});
      } else if ( max_i == intrvl[1] ){
	curves.push_back( { {series[max_i], 0}, intrvl[1]} );
	stack.push_back( {intrvl[0],intrvl[1]-1});
      } else {
	DoublePair s1_curve = pla::regression( series.data() + intrvl[0], series.data() + max_i - 1);
	DoublePair s2_curve = pla::regression( series.data() + max_i + 1, series.data() + intrvl[1]);

	bool add_break_to_s1 = dist_from_line(s1_curve, double(max_i-intrvl[0]), series[max_i]) 
				< dist_from_line(s2_curve, -1.0, series[max_i]);
	if (add_break_to_s1) {
	  stack.push_back( {intrvl[0], max_i});
	  stack.push_back( {max_i + 1,intrvl[1]});
	} else {
	  stack.push_back( {intrvl[0],max_i - 1});
	  stack.push_back( {max_i,intrvl[1]});
	}
      }
    }

  }
  std::sort(curves.begin(), curves.end(), [](const auto& tp1, const auto& tp2) { return std::get<1>(tp1) < std::get<1>(tp2); });
  return curves;
}

vector<tuple<DoublePair, unsigned int>> dac_curve_fitting::dac_linear_early_cutoff( const vector<double>& series, double epsilon, unsigned int num_seg)
{
  if (epsilon < 0) return {};
  if (series.size() == 0) return {};

  vector<array<unsigned int, 2>> queue = {};
  queue.push_back( {0, (unsigned int) series.size() - 1} );
  vector<tuple<DoublePair, unsigned int>> curves;

  auto dist_from_line = [](DoublePair line, double x, double y) { return std::abs( y - line[0] - line[1]*x); };

  while (queue.size() > 0 && curves.size() < num_seg) {
    array<unsigned int, 2> intrvl = queue.back();
    queue.pop_back();
    DoublePair regressed = pla::regression( series.data() + intrvl[0], series.data() + intrvl[1] );
    
    unsigned int max_i = intrvl[0];
    double max_dist = 0;
    for (int i=intrvl[0]; i<=intrvl[1]; ++i) {
      if ( dist_from_line(regressed, i-intrvl[0], series[i]) >= dist_from_line( regressed, max_i-intrvl[0], series[max_i])) {
	max_i = i;
	max_dist = dist_from_line(regressed, i-intrvl[0], series[i]);
      }
    }

    if ( max_dist <= epsilon) {
      curves.push_back( {regressed, intrvl[1] } );
    } else {
      if ( max_i == intrvl[0]) {
	curves.push_back( { {series[max_i], 0 }, intrvl[0]} );
	queue.push_back( {intrvl[0]+1,intrvl[1]});
      } else if ( max_i == intrvl[1] ){
	curves.push_back( { {series[max_i], 0}, intrvl[1]} );
	queue.push_back( {intrvl[0],intrvl[1]-1});
      } else {
	DoublePair s1_curve = pla::regression( series.data() + intrvl[0], series.data() + max_i - 1);
	DoublePair s2_curve = pla::regression( series.data() + max_i + 1, series.data() + intrvl[1]);

	bool add_break_to_s1 = dist_from_line(s1_curve, double(max_i-intrvl[0]), series[max_i]) < dist_from_line(s2_curve, -1.0, series[max_i]);
	if (add_break_to_s1) {
	  queue.push_back( {intrvl[0], max_i});
	  queue.push_back( {max_i + 1,intrvl[1]});
	} else {
	  queue.push_back( {intrvl[0],max_i - 1});
	  queue.push_back( {max_i,intrvl[1]});
	}
      }
    }

  }
  std::sort(curves.begin(), curves.end(), [](const auto& tp1, const auto& tp2) { return std::get<1>(tp1) < std::get<1>(tp2); });
  return curves;
}
