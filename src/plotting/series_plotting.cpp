#include "series_plotting.h"

#include "gnuplot-iostream.h"
#include "mse.h"
#include "plotting_constants.h"
#include <boost/tuple/tuple.hpp>

#include "ucr_parsing.h"
#include "z_norm.h"

using std::vector;

void plot_series(Series& s1, std::string file_path)
{
  vector<double> x,y;
  int i=0;
  for (const auto& f : s1.series) {
    x.emplace_back( (double) i);
    y.emplace_back( (double) f);
    ++i;
  }
  Gnuplot gp;
  gp << "set term pdfcairo enhanced color dashed font 'Verdana, 14' rounded size 32cm, 19.2cm\n";
  gp << "set output '"<< file_path << s1.name <<".pdf' \n";
  gp << "set xlabel 'time'\n";
  gp << "set ylabel 'series value'\n";
  gp << "set title 'Plot of series " << s1.name << "'\n";
  gp << "plot '-' with linespoints lt rgb 'blue' lw 1.2 pt 5 ps 0.5 title '" <<  s1.name << "'\n";
  gp.send1d( boost::make_tuple( x, y ) );

  std::string command = "firefox '" + file_path + s1.name + ".pdf'";
  system(command.c_str()); 
}

void plot_many_series(vector<Series>& vs)
{
  if (vs.size() ==0) return;

  vector<vector<double>> vx,vy;
  int i=0;
  for (int s_i=0; s_i < vs.size(); ++s_i) {
    i=0;
    vx.push_back({});
    vy.push_back({});
    for (const auto& f : vs[s_i].series ) {
      vx[s_i].emplace_back( (double) i);
      vy[s_i].emplace_back( (double) f);
      i++;
    }
  }
  Gnuplot gp;
  using namespace gp_constants;
  gp << "set terminal " << GNUPLOT_TERMINAL << "\n";
  gp << "set term " << GNUPLOT_TERMINAL <<" size " << GNUPLOT_SIZE_X << " " << GNUPLOT_SIZE_Y << "\n";
  gp << "set xlabel 'time'\n";
  gp << "set ylabel 'series val'\n";
  gp << "set title 'Plot of multiple time series \n";
  gp << "plot ";
  for (int i=0; i<vs.size()-1; ++i) {
    gp << "'-' with lines title '" <<  vs[i].name << "', ";
  }
  gp << "'-' with lines title '" <<  vs.back().name << "'\n";
  for (int i=0; i<vs.size(); ++i) {
    gp.send1d( boost::make_tuple( vx[i], vy[i] ) );
  }
}

int min(int x, int y) { if (x<y) return x; return y;}

void plot_series_diff(Series& s1, Series& s2, std::string file_path)
{
  vector<double> x1,y1,x2,y2,err1,err2;
  int i=0;
  for (int i=0; i<min(s1.series.size(), s2.series.size()); ++i ) {
    x1.emplace_back( (double) i);
    y1.emplace_back( (double) s1.series[i]);
    x2.emplace_back( (double) i);
    y2.emplace_back( (double) s2.series[i]);
    err1.emplace_back(s2.series[i]);
    err2.emplace_back(s1.series[i]);
    ++i;
  }
  Gnuplot gp;
  gp << "set term pdfcairo enhanced color dashed font 'Verdana, 14' rounded size 32cm, 19.2cm\n";
  gp << "set output '"<< file_path << s1.name <<" and " << s2.name << ".pdf' \n";
  gp << "set xlabel 'time'\n";
  gp << "set ylabel 'series values'\n";
  gp << "set title 'Plot of series " << s1.name << " and series "<< s2.name <<"'\n";
  gp << "plot '-' dt 3 lw 2 with yerrorlines title 'difference'"
     << ", '-' with linespoints lt rgb 'blue' lw 1.2 pt 5 ps 0.5 title '" << s1.name << "'"
     << ", '-' with linespoints lt rgb 'red' lw 1.2 pt 7 ps 0.5 title '" << s2.name << "'\n";
  gp.send1d( boost::make_tuple( x1, y1, err1, err2 ) );
  gp.send1d( boost::make_tuple( x1, y1 ) );
  gp.send1d( boost::make_tuple( x2, y2 ) );

  std::string command = "firefox '" + file_path + s1.name + " and " + s2.name + ".pdf'";
  system(command.c_str()); 
}

void plot_lines(vector<Line> lines, PlotLabels pl)
{

  if (lines.size() ==0) return;

  Gnuplot gp;
  using namespace gp_constants;
  // gp << "set terminal " << GNUPLOT_TERMINAL << "\n";
  // gp << "set term " << GNUPLOT_TERMINAL <<" size " << GNUPLOT_SIZE_X << " " << GNUPLOT_SIZE_Y << "\n";
  gp << "set term png size 1280 640 background 'white' enhanced font size 20 \n";
  gp << "set output 'img/drt_comparisons/"<<pl.title<<".png' \n";
  gp << "set xlabel '" << pl.xlabel << "'\n";
  gp << "set ylabel '" << pl.ylabel << "'\n";
  gp << "set title '" << pl.title << "'\n";
  gp << "plot ";
  for (int i=0; i<lines.size()-1; ++i) {
    gp << "'-' with lines title '" <<  lines[i].name << "', ";
  }
  gp << "'-' with lines title '" <<  lines.back().name << "'\n";
  for (int i=0; i<lines.size(); ++i) {
    gp.send1d( boost::make_tuple( lines[i].x, lines[i].y ) );
  }
}

void plot_ucr_drt_mse(const vector<std::string>& datasets
		      , std::string datasets_loc
		      , unsigned int ds_start
		      , unsigned int ds_end
		      , unsigned int max_num_params
		      , std::function<std::vector<double> (const std::vector<double>&, unsigned int num_params)> drt)
{
  if (max_num_params < 5) return;

  vector<double> x;
  for (int j=4; j<=max_num_params; j+=100) x.push_back(j);

  vector<Line> lines;
  vector<vector<double>> datasets_mse(ds_end-ds_start+1);
  for (int i=ds_start; i<=ds_end; ++i) {
    vector<double> dataset_i = ucr_parsing::parse_ucr_dataset(datasets[i], datasets_loc, ucr_parsing::DatasetType::TRAIN);
    z_norm::z_normalise(dataset_i);
    for (int j=4; j<=max_num_params; j+=100) {
      datasets_mse[i-ds_start].push_back( mse::mse_between_seq(dataset_i, drt(dataset_i, j)) );
      std::cout << mse::mse_between_seq(dataset_i, drt(dataset_i,j)) << std::endl;
    }
    lines.push_back( { x, datasets_mse[i-ds_start], "dataset " + std::to_string(i) +": "+datasets[i] } );
  }
  PlotLabels pl = {"comparison of MSE for DRT on "+ std::to_string(ds_end-ds_start+1) +" datasets", "number of parameters", "MSE error"};
  plot_lines( lines, pl);
}

void plot_ucr_drt_maxdev(const vector<std::string>& datasets
		      , std::string datasets_loc
		      , unsigned int ds_start
		      , unsigned int ds_end
		      , unsigned int max_num_params
		      , std::function<std::vector<double> (const std::vector<double>&, unsigned int num_params)> drt)
{
  if (max_num_params < 5) return;

  vector<double> x;
  for (int j=4; j<=max_num_params; j+=100) x.push_back(j);

  vector<Line> lines;
  vector<vector<double>> datasets_mse(ds_end-ds_start+1);
  for (int i=ds_start; i<=ds_end; ++i) {
    vector<double> dataset_i = ucr_parsing::parse_ucr_dataset(datasets[i], datasets_loc, ucr_parsing::DatasetType::TRAIN);
    z_norm::z_normalise(dataset_i);
    for (int j=4; j<=max_num_params; j+=100) {
      datasets_mse[i-ds_start].push_back( mse::maxdev_between_seq(dataset_i, drt(dataset_i, j)) );
    }
    lines.push_back( { x, datasets_mse[i-ds_start], "dataset " + std::to_string(i) +": "+datasets[i] } );
  }
  PlotLabels pl = {"comparison of maximum deviation for DRT on "+ std::to_string(ds_end-ds_start+1) +" datasets", "number of parameters", "MaxDev error"};
  plot_lines( lines, pl);
}
