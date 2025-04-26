#include "plot_dimreduct_pla.h"

#include <functional>
#include <math.h>
#include <tuple>

#include "gnuplot-iostream.h"
#include <boost/tuple/tuple.hpp>

/**
 * @file plot_dimreduct_pla.cpp
 * @brief Provides implementations for plot_dimreduct_pla.h .
 */

using std::vector;
using std::tuple;


void plot_pla::plot_pla(const vector<double>& series, std::string data_name, unsigned int num_params, PlotDetails p)
{
  if (num_params <= 0) return;

  vector<DoublePair> pla_subseq = pla::pla(series, num_params);

  unsigned int interval_size = (series.size() / pla_subseq.size()) + (series.size() % pla_subseq.size() != 0);

  vector<double> x1, x2, err1, err2, y1, y2;
  for (int i=0; i<series.size(); ++i){
    x1.emplace_back(i);
    x2.emplace_back(i);

    y1.emplace_back( series[i] );
    y2.emplace_back( pla_subseq[i / interval_size][0] + pla_subseq[i / interval_size][1] * (i % interval_size));

    err1.emplace_back( series[i] );
    err2.emplace_back( pla_subseq[i / interval_size][0] + pla_subseq[i / interval_size][1] * (i % interval_size));

    if ( (i+1) % interval_size == 0) {
      x2.emplace_back( (double)i + 0.50);
      y2.emplace_back( NAN );
    }
  }

  Gnuplot gp;
  plot_setup::setup_gnuplot(gp, p);
  //gp << "plot '-' dt 3 with yerrorlines title 'errors', '-' with lines lt rgb 'blue' lw 1.2 title '"<<data_name<<"', '-' with linespoints lt rgb 'red' lw 1.2 title 'PLA'\n";
  gp << "plot '-' with lines lt rgb 'blue' lw 4 title '"<<data_name<<"', '-' with lines lt rgb 'red' lw 4 title 'PLA'\n";
  //gp.send1d( boost::make_tuple( x1, y1, err1, err2) );
  gp.send1d( boost::make_tuple( x1, y1 ) );
  gp.send1d( boost::make_tuple( x2, y2 ) );
  plot_setup::open_pdf(p);
}

void plot_pla::plot_any_apla(const vector<double>& series, std::string data_name
			  , std::function<const Seqddt(const Seqd&)> apla_func
			  , std::string drt_name
			  , PlotDetails p)
{
  vector<tuple<DoublePair, unsigned int>> apla_series = apla_func(series);

  vector<double> x1, x2, err1, err2, y1, y2;
  int i=0;
  int j=0;
  for (const auto& [dp,end_i] : apla_series) {
    j=0;
    while ( i <= end_i) {
      x1.emplace_back(i);
      x2.emplace_back(i);

      y1.emplace_back( series[i] );
      y2.emplace_back( dp[0] + dp[1]*j );

      err1.emplace_back( series[i] );
      err2.emplace_back( dp[0] + dp[1]*j );

      if ( (i+1) > end_i) {
	x2.emplace_back( (double)i + 0.50);
	y2.emplace_back( NAN );
      }
      i++;
      j++;
    }
  }

  Gnuplot gp;
  plot_setup::setup_gnuplot(gp, p);
  gp << "plot '-' dt 3 with yerrorlines title 'Errors', '-' with linespoints lt rgb 'blue' lw 1.2 pt 5 ps 0.5 title '"<<data_name<<"', '-' with linespoints lt rgb 'red' lw 1.2 pt 7 ps 0.5 title '"<< drt_name <<"'\n";
  //gp << "plot '-' with lines lt rgb 'blue' lw 4 title '"<<data_name<<"', '-' with lines lt rgb 'red' lw 4 title '"<< drt_name <<"'\n";
  gp.send1d( boost::make_tuple( x1, y1, err1, err2) );
  gp.send1d( boost::make_tuple( x1, y1 ) );
  gp.send1d( boost::make_tuple( x2, y2 ) );
  plot_setup::open_pdf(p);
}
