#include "dac_curve_fitting.h"

#include <numeric>
#include <algorithm>
#include <cmath>
#include <array>

#include <iostream>
using std::cout, std::endl;

using std::vector;
using std::tuple;
using std::array;

// uses an IEEE float trick to return NaN, won't work if compiler doesn't support IEEE
float mean(const float* const f1, const float* const f2)
{
  if (f2 < f1) return 0.0 / 0.0; // return NaN
  int size = f2 - f1 + 1;
  return std::accumulate(f1, f2+1, 0.0) / (float) size;
}

vector<tuple<float, unsigned int>> dac_curve_fitting::dac_constant( const vector<float>& series, float epsilon)
{
  if (epsilon < 0) return {};
  if (series.size() == 0) return {};

  vector<array<int, 2>> queue = {};
  queue.push_back( {0, (int) series.size() - 1} );
  vector<tuple<float, unsigned int>> curves;

  while (queue.size() > 0) {
    array<int, 2> intrvl = queue.back();
    queue.pop_back();
    float s_mean = mean( series.data() + intrvl[0], series.data() + intrvl[1] );
    
    const auto max_dev_el = std::max_element(series.begin() + intrvl[0], series.begin() + intrvl[1] + 1,
	[&s_mean](float i, float j){ return std::abs(i-s_mean) < std::abs(j-s_mean); }
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

	float s1_mean = mean( series.data() + intrvl[0], series.data() + intrvl[0] + break_point - 1);
	float s2_mean = mean( series.data() + intrvl[0] + break_point + 1, series.data() + intrvl[1]);

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
