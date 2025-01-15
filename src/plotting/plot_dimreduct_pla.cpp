#include "pla.h"
#include "dac_curve_fitting.h"

#include <string>
#include <algorithm>
#include <vector>
#include <tuple>
#include <functional>

#include <iostream>

#include <math.h>

#include "gnuplot-iostream.h"
#include "plotting_constants.h"
#include <boost/tuple/tuple.hpp>

#include "series_plotting.h"
#include "ucr_parsing.h"
#include "z_norm.h"

using std::vector;
using std::tuple;


void plot_pla_subseq(const vector<double>& series, std::string data_name, unsigned int num_params, unsigned int start_pos, unsigned int end_pos)
{
  if (start_pos >= end_pos) return;
  if (num_params <= 0) return;

  vector<double> subseq(end_pos - start_pos + 1);
  std::copy(series.cbegin() + start_pos, series.cbegin() + end_pos + 1, subseq.begin());

  vector<DoublePair> pla_subseq = pla::chunk_regression(subseq, num_params);
  std::cout << " PLA uses " << pla_subseq.size() << " intervals " << std::endl;

  unsigned int interval_size = (subseq.size() / pla_subseq.size()) + (subseq.size() % pla_subseq.size() != 0);

  vector<double> x1, x2, err1, err2, y1, y2;
  for (int i=0; i<subseq.size(); ++i){
    x1.emplace_back(i);
    x2.emplace_back(i);

    y1.emplace_back( subseq[i] );
    y2.emplace_back( pla_subseq[i / interval_size][0] + pla_subseq[i / interval_size][1] * (i % interval_size));

    err1.emplace_back( subseq[i] );
    err2.emplace_back( pla_subseq[i / interval_size][0] + pla_subseq[i / interval_size][1] * (i % interval_size));

    if ( (i+1) % interval_size == 0) {
      x2.emplace_back( (double)i + 0.50);
      y2.emplace_back( NAN );
    }
  }

  Gnuplot gp;
  using namespace gp_constants;
  // gp << "set terminal " << GNUPLOT_TERMINAL << "\n";
  // gp << "set term " << GNUPLOT_TERMINAL <<" size " << GNUPLOT_SIZE_X << " " << GNUPLOT_SIZE_Y << "\n";
  gp << "set term png size 1280 640 background 'white' enhanced font size 20 \n";
  gp << "set output 'img/drt_comparisons/"<<data_name<<"_pla_subseq_"<<num_params<<"_params.png' \n";
  gp << "set xlabel 'time'\n";
  gp << "set ylabel 'series val'\n";
  gp << "set title 'Plot of "<<data_name<<" against PLA approximation'\n";
  gp << "plot '-' dt 3 with yerrorlines title 'errors', '-' with lines lt rgb 'blue' lw 1.2 title '"<<data_name<<"', '-' with linespoints lt rgb 'red' lw 1.2 title 'PLA'\n";
  gp.send1d( boost::make_tuple( x1, y1, err1, err2) );
  gp.send1d( boost::make_tuple( x1, y1 ) );
  gp.send1d( boost::make_tuple( x2, y2 ) );
}

void plot_pla_mse(const vector<std::string>& datasets, std::string datasets_loc, int start, int end, int max_int_size)
{
  vector<double> x;
  for (int j=2; j<=max_int_size; ++j) x.push_back(j);

  vector<Line> lines;
  vector<vector<double>> datasets_mse(end-start+1);
  for (int i=start; i<=end; ++i) {
    vector<double> dataset_i = ucr_parsing::parse_ucr_dataset(datasets[i], datasets_loc, ucr_parsing::DatasetType::TEST);
    z_norm::z_normalise(dataset_i);
    for (int j=2; j<=max_int_size; ++j) {
      datasets_mse[i-start].push_back( pla::pla_mse(dataset_i, j) );
    }
    lines.push_back( { x, datasets_mse[i-start], "dataset " + std::to_string(i) +": "+datasets[i] } );
  }
  PlotLabels pl = {"comparison of pla on "+ std::to_string(end-start+1) +" datasets", "interval size", "MSE error"};
  plot_lines( lines, pl);
}

void plot_dac_apla_subseq(const vector<double>& series, std::string data_name, double epsilon, unsigned int start_pos, unsigned int end_pos)
{
  using std::tuple;
  if (start_pos >= end_pos) return;
  if (epsilon < 0) return;

  vector<double> subseq(end_pos - start_pos + 1);
  std::copy(series.cbegin() + start_pos, series.cbegin() + end_pos + 1, subseq.begin());

  vector<tuple<DoublePair, unsigned int>> apla_subseq = dac_curve_fitting::dac_linear(subseq, epsilon);
  std::cout << " APLA uses " << apla_subseq.size() << " intervals " << std::endl;

  vector<double> x1, x2, err1, err2, y1, y2;
  int i=0;
  int j=0;
  for (const auto& tp : apla_subseq) {
    j=0;
    while ( i <= std::get<1>(tp)) {
      x1.emplace_back(i+start_pos);
      x2.emplace_back(i+start_pos);

      y1.emplace_back( subseq[i] );
      y2.emplace_back( std::get<0>(tp)[0] + std::get<0>(tp)[1]*j );

      err1.emplace_back( subseq[i] );
      err2.emplace_back( std::get<0>(tp)[0] + std::get<0>(tp)[1]*j );

      if ( (i+1) > std::get<1>(tp)) {
	x2.emplace_back( (double)i + 0.50);
	y2.emplace_back( NAN );
      }
      i++;
      j++;
    }
  }

  Gnuplot gp;
  using namespace gp_constants;
  //gp << "set terminal " << GNUPLOT_TERMINAL << "\n";
  //gp << "set term " << GNUPLOT_TERMINAL <<" size " << GNUPLOT_SIZE_X << " " << GNUPLOT_SIZE_Y << "\n";
  gp << "set term png size 1280 640 background 'white' enhanced font size 20 \n";
  gp << "set output 'img/drt_comparisons/"<<data_name<<"_rec_apla_subseq_"<<epsilon<<"_deviation.png' \n";
  gp << "set xlabel 'time'\n";
  gp << "set ylabel 'series val'\n";
  gp << "set title 'Plot of "<<data_name<<" against APLA approximation'\n";
  gp << "plot '-' dt 3 with yerrorlines title 'errors', '-' with lines lt rgb 'blue' lw 1.2 title '"<<data_name<<"', '-' with linespoints lt rgb 'red' lw 1.2 title 'APLA'\n";
  gp.send1d( boost::make_tuple( x1, y1, err1, err2) );
  gp.send1d( boost::make_tuple( x1, y1 ) );
  gp.send1d( boost::make_tuple( x2, y2 ) );
}

void plot_any_apla_subseq(const vector<double>& series, std::string data_name
			  , std::function<const vector<tuple<DoublePair, unsigned int>>(const vector<double>&)> apla_func
			  , unsigned int start_pos, unsigned int end_pos)
{
  using std::tuple;
  if (start_pos >= end_pos) return;

  vector<double> subseq(end_pos - start_pos + 1);
  std::copy(series.cbegin() + start_pos, series.cbegin() + end_pos + 1, subseq.begin());

  vector<tuple<DoublePair, unsigned int>> apla_subseq = apla_func(subseq);
  std::cout << " APLA uses " << apla_subseq.size() << " intervals " << std::endl;

  vector<double> x1, x2, err1, err2, y1, y2;
  int i=0;
  int j=0;
  for (const auto& tp : apla_subseq) {
    j=0;
    while ( i <= std::get<1>(tp)) {
      x1.emplace_back(i+start_pos);
      x2.emplace_back(i+start_pos);

      y1.emplace_back( subseq[i] );
      y2.emplace_back( std::get<0>(tp)[0] + std::get<0>(tp)[1]*j );

      err1.emplace_back( subseq[i] );
      err2.emplace_back( std::get<0>(tp)[0] + std::get<0>(tp)[1]*j );

      if ( (i+1) > std::get<1>(tp)) {
	x2.emplace_back( (double)i + 0.50);
	y2.emplace_back( NAN );
      }
      i++;
      j++;
    }
  }

  Gnuplot gp;
  using namespace gp_constants;
  gp << "set terminal " << GNUPLOT_TERMINAL << "\n";
  gp << "set term " << GNUPLOT_TERMINAL <<" size " << GNUPLOT_SIZE_X << " " << GNUPLOT_SIZE_Y << "\n";
  // gp << "set term png size 1280 640 background 'white' enhanced font size 20 \n";
  // gp << "set output 'img/drt_comparisons/"<<data_name<<"_apla_subseq.png' \n";
  gp << "set xlabel 'time'\n";
  gp << "set ylabel 'series val'\n";
  gp << "set title 'Plot of "<<data_name<<" against APLA approximation'\n";
  gp << "plot '-' dt 3 with yerrorlines title 'errors', '-' with lines lt rgb 'blue' lw 1.2 title '"<<data_name<<"', '-' with linespoints lt rgb 'red' lw 1.2 title 'APLA'\n";
  gp.send1d( boost::make_tuple( x1, y1, err1, err2) );
  gp.send1d( boost::make_tuple( x1, y1 ) );
  gp.send1d( boost::make_tuple( x2, y2 ) );
}
