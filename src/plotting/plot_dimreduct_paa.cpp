#include "paa.h"
#include "dac_curve_fitting.h"

#include <string>
#include <algorithm>
#include <vector>

#include <math.h>

#include <iostream>

#include "gnuplot-iostream.h"
#include "plotting_constants.h"
#include <boost/tuple/tuple.hpp>

#include "series_plotting.h"
#include "ucr_parsing.h"
#include "z_norm.h"

using std::vector;


void plot_paa_subseq_series_comparison(const vector<double>& series, std::string data_name, unsigned int interval_size, unsigned int start_pos, unsigned int end_pos)
{
  if (start_pos >= end_pos) return;
  if (interval_size <= 0) return;

  vector<double> subseq(end_pos - start_pos + 1);
  std::copy(series.cbegin() + start_pos, series.cbegin() + end_pos + 1, subseq.begin());

  vector<double> paa_subseq = paa::paa(subseq, interval_size);
  std::cout << "PAA uses " << paa_subseq.size() << " intervals" << std::endl;

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
  using namespace gp_constants;
  gp << "set terminal " << GNUPLOT_TERMINAL << "\n";
  gp << "set term " << GNUPLOT_TERMINAL <<" size " << GNUPLOT_SIZE_X << " " << GNUPLOT_SIZE_Y << "\n";
  gp << "set xlabel 'time'\n";
  gp << "set ylabel 'series val'\n";
  gp << "set title 'Plot of series against PAA approximation'\n";
  gp << "plot '-' dt 3 with yerrorlines title 'errors', '-' with lines title 'Series', '-' with linespoints title 'PAA'\n";
  gp.send1d( boost::make_tuple( x1, y1, err1, err2) );
  gp.send1d( boost::make_tuple( x1, y1 ) );
  gp.send1d( boost::make_tuple( x2, y2 ) );
}

void plot_paa_mse(const vector<std::string>& datasets, std::string datasets_loc, int start, int end, int max_int_size)
{
  vector<double> x;
  for (int j=1; j<=max_int_size; ++j) x.push_back(j);

  vector<Line> lines;
  vector<vector<double>> datasets_mse(end-start+1);
  for (int i=start; i<=end; ++i) {
    vector<double> dataset_i = ucr_parsing::parse_ucr_dataset(datasets[i], datasets_loc, ucr_parsing::DatasetType::TEST);
    z_norm::z_normalise(dataset_i);
    for (int j=1; j<=max_int_size; ++j) {
      datasets_mse[i-start].push_back( paa::paa_mse(dataset_i, j) );
    }
    lines.push_back( { x, datasets_mse[i-start], "dataset " + std::to_string(i) +": "+datasets[i] } );
  }
  PlotLabels pl = {"comparison of paa on "+ std::to_string(end-start+1) +" datasets", "interval size", "MSE error"};
  plot_lines( lines, pl);
}

void plot_apaa_subseq_series_comparison(const vector<double>& series, std::string data_name, double epsilon, unsigned int start_pos, unsigned int end_pos)
{
  using std::tuple;
  if (start_pos >= end_pos) return;
  if (epsilon < 0) return;

  vector<double> subseq(end_pos - start_pos + 1);
  std::copy(series.cbegin() + start_pos, series.cbegin() + end_pos + 1, subseq.begin());

  vector<tuple<double, unsigned int>> apaa_subseq = dac_curve_fitting::dac_constant(subseq, epsilon);
  std::cout << "APAA uses " << apaa_subseq.size() << " intervals" << std::endl;

  vector<double> x1, x2, err1, err2, y1, y2;
  int i=0;
  for (const auto& tp : apaa_subseq) {
    while ( i <= std::get<1>(tp)) {
      x1.emplace_back(i+start_pos);
      x2.emplace_back(i+start_pos);

      y1.emplace_back( subseq[i] );
      y2.emplace_back( std::get<0>(tp) );

      err1.emplace_back( subseq[i] );
      err2.emplace_back( std::get<0>(tp));

      if ( (i+1) > std::get<1>(tp)) {
	x1.emplace_back( (double)i + 0.50);
	y1.emplace_back( NAN );
	err1.emplace_back( subseq[i] );
	err2.emplace_back( subseq[i] );

	x2.emplace_back( (double)i + 0.50);
	y2.emplace_back( NAN );
      }
      i++;
    }
  }

  Gnuplot gp;
  using namespace gp_constants;
  gp << "set terminal " << GNUPLOT_TERMINAL << "\n";
  gp << "set term " << GNUPLOT_TERMINAL <<" size " << GNUPLOT_SIZE_X << " " << GNUPLOT_SIZE_Y << "\n";
  gp << "set xlabel 'time'\n";
  gp << "set ylabel 'series val'\n";
  gp << "set title 'Plot of "<<data_name<<" against APAA approximation'\n";
  gp << "plot '-' dt 3 with yerrorlines title 'errors', '-' with lines title '"<<data_name<<"', '-' with linespoints title 'APAA'\n";
  gp.send1d( boost::make_tuple( x1, y1, err1, err2) );
  gp.send1d( boost::make_tuple( x1, y1 ) );
  gp.send1d( boost::make_tuple( x2, y2 ) );
}
