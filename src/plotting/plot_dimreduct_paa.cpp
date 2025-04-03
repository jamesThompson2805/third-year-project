#include "plot_dimreduct_paa.h"

/**
 * @file plot_dimreduct_paa.cpp
 * @brief Source file containing implementations of functions defined in plot_dimreduct_paa.h .
 */

#include "paa.h"

#include <math.h>


#include "gnuplot-iostream.h"
#include "plot_types.h"
#include <boost/tuple/tuple.hpp>

using std::vector;
using std::tuple;

/**
 * @brief plot_paa plots a series with its PAA approximation render on top.
 * @param series Is the time series you want plotted in entirity.
 * @param data_name Is the name of the data passed as series.
 * @param num_params Is the number of parameters allocated to PAA, this is identical to the number of segments PAA will use.
 * @param p Is the PlotDetails for the rendered plot.
 */
void plot_paa::plot_paa(const vector<double>& series, std::string data_name, unsigned int num_params, PlotDetails p)
{
  if (num_params <= 0) return;

  vector<double> paa = paa::paa(series, num_params);

  unsigned int interval_size = (series.size() / paa.size()) + (series.size() % paa.size() != 0);

  vector<double> x1, x2, err1, err2, y1, y2;
  for (int i=0; i<series.size(); ++i){
    x1.emplace_back(i);
    x2.emplace_back(i);

    y1.emplace_back( series[i] );
    y2.emplace_back( paa[i / interval_size] );

    err1.emplace_back( series[i] );
    err2.emplace_back( paa[ i / interval_size]);

    if ( (i+1) % interval_size == 0) {
      x2.emplace_back( (double)i + 0.50);
      y2.emplace_back( NAN );
    }
  }

  Gnuplot gp;
  plot_setup::setup_gnuplot(gp, p);
  //gp << "plot '-' dt 3 with yerrorlines title 'errors', '-' with lines lt rgb 'blue' lw 1.2 title '"<<data_name<<"', '-' with linespoints lt rgb 'red' lw 1.2 title 'PAA'\n";
  gp << "plot '-' with lines lt rgb 'blue' lw 4 title '"<<data_name<<"', '-' with lines lt rgb 'red' lw 4 title 'PAA'\n";
  //gp.send1d( boost::make_tuple( x1, y1, err1, err2) );
  gp.send1d( boost::make_tuple( x1, y1 ) );
  gp.send1d( boost::make_tuple( x2, y2 ) );
  plot_setup::open_pdf(p);
}

/**
 * @brief plot_any_apaa plots a series over time with its approximation by an adaptive PAA implementation on top.
 * @param series Is the time series you want plotted in entirity.
 * @param data_name Is the name of the data passed as series.
 * @param to_apaa Is a function that takes the series you pass in and returns an adaptive PAA representation of it. The form recognised is a tuple of double and unsigned int where the double is the value upon the region defined by previous index to the current index (the unsigned integer).
 * @param drt_name Is the name of the dimension reduction technique passed eg. APCA.
 * @param p Is the PlotDetails for the rendered plot.
 */
void plot_paa::plot_any_apaa(const vector<double>& series
			  , std::string data_name
			  , std::function< vector<tuple<double, unsigned int>>(const vector<double>&) > to_apaa
			  , std::string drt_name
			  , PlotDetails p)
{

  vector<tuple<double, unsigned int>> apaa = to_apaa(series);

  vector<double> x1, x2, err1, err2, y1, y2;
  int i=0;
  for (const auto& tp : apaa) {
    while ( i <= std::get<1>(tp)) {
      x1.emplace_back(i);
      x2.emplace_back(i);

      y1.emplace_back( series[i] );
      y2.emplace_back( std::get<0>(tp) );

      err1.emplace_back( series[i] );
      err2.emplace_back( std::get<0>(tp));

      if ( (i+1) > std::get<1>(tp)) {
	x2.emplace_back( (double)i + 0.50);
	y2.emplace_back( NAN );
      }
      i++;
    }
  }

  Gnuplot gp;
  plot_setup::setup_gnuplot(gp, p);
  //gp << "plot '-' dt 3 with yerrorlines title 'Errors', '-' with linespoints lt rgb 'blue' lw 1.2 pt 5 ps 0.5 title '"<<data_name<<"', '-' with linespoints lt rgb 'red' lw 1.2 pt 7 ps 0.5 title '"<< drt_name <<"'\n";
  gp << "plot '-' with lines lt rgb 'blue' lw 4 title '"<<data_name<<"', '-' with lines lt rgb 'red' lw 4 title '"<< drt_name <<"'\n";
  //gp.send1d( boost::make_tuple( x1, y1, err1, err2) );
  gp.send1d( boost::make_tuple( x1, y1 ) );
  gp.send1d( boost::make_tuple( x2, y2 ) );

  plot_setup::open_pdf(p);
}
