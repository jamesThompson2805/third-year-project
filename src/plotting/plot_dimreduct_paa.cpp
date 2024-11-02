#include "paa.h"

#include <string>
#include <algorithm>
#include <vector>

#include <math.h>

#include "gnuplot-iostream.h"
#include <boost/tuple/tuple.hpp>

using std::vector;


void plot_subseq_series_comparison(const vector<float>& series, std::string data_name, unsigned int interval_size, unsigned int start_pos, unsigned int end_pos)
{
  if (start_pos >= end_pos) return;
  if (interval_size <= 0) return;

  vector<float> subseq(end_pos - start_pos + 1);
  std::copy(series.cbegin() + start_pos, series.cbegin() + end_pos + 1, subseq.begin());

  vector<float> paa_subseq = paa::paa(subseq, interval_size);

  vector<double> x1, x2, err1, err2, y1, y2;
  for (int i=0; i<subseq.size(); ++i){
    x1.emplace_back(i);
    x2.emplace_back(i);

    y1.emplace_back( subseq[i] );
    y2.emplace_back( paa_subseq[i / interval_size] );

    err1.emplace_back( subseq[i] );
    err2.emplace_back( paa_subseq[ i / interval_size]);

    if ( (i+1) % interval_size == 0) {
      x1.emplace_back( (double)i + 0.50);
      y1.emplace_back( NAN );
      err1.emplace_back( subseq[i] );
      err2.emplace_back( subseq[i] );

      x2.emplace_back( (double)i + 0.50);
      y2.emplace_back( NAN );
    }
  }

  Gnuplot gp;
  gp << "set terminal x11\n";
  gp << "set term x11 size 1280 800 \n";
  gp << "set xlabel 'time'\n";
  gp << "set ylabel 'series val'\n";
  gp << "set title 'Plot of series against PAA approximation'\n";
  gp << "plot '-' dt 3 with yerrorlines title 'errors', '-' with lines title 'Series', '-' with linespoints title 'PAA'\n";
  gp.send1d( boost::make_tuple( x1, y1, err1, err2) );
  gp.send1d( boost::make_tuple( x1, y1 ) );
  gp.send1d( boost::make_tuple( x2, y2 ) );
}
