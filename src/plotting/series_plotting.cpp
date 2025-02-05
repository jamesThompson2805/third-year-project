#include "series_plotting.h"

#include "gnuplot-iostream.h"
#include "mse.h"
#include "plot_types.h"
#include "plotting_constants.h"
#include <boost/tuple/tuple.hpp>

#include "ucr_parsing.h"
#include "z_norm.h"

using std::vector;

void plot::plot_series(Series& s1, PlotDetails p, Backend b)
{
  vector<double> x,y;
  int i=0;
  for (const auto& f : s1.series) {
    x.emplace_back( (double) i);
    y.emplace_back( (double) f);
    ++i;
  }
  Gnuplot gp;
  plot_setup::setup_gnuplot(gp, b, p);
  gp << "plot '-' with linespoints lt rgb 'blue' lw 1.2 pt 5 ps 0.5 title '" <<  s1.name << "'\n";
  gp.send1d( boost::make_tuple( x, y ) );

  if (b == PDF)
    plot_setup::open_pdf(p);
}

void plot::plot_many_series(vector<Series>& vs, PlotDetails p, Backend b)
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
  plot_setup::setup_gnuplot(gp, b, p);
  gp << "plot ";
  for (int i=0; i<vs.size()-1; ++i) {
    gp << "'-' with lines title '" <<  vs[i].name << "', ";
  }
  gp << "'-' with lines title '" <<  vs.back().name << "'\n";
  for (int i=0; i<vs.size(); ++i) {
    gp.send1d( boost::make_tuple( vx[i], vy[i] ) );
  }

  if (b == PDF)
    plot_setup::open_pdf(p);
}

int min(int x, int y) { if (x<y) return x; return y;}

void plot::plot_series_diff(Series& s1, Series& s2, PlotDetails p, Backend b)
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
  plot_setup::setup_gnuplot(gp, b, p);
  gp << "plot '-' dt 3 lw 2 with yerrorlines title 'difference'"
     << ", '-' with linespoints lt rgb 'blue' lw 1.2 pt 5 ps 0.5 title '" << s1.name << "'"
     << ", '-' with linespoints lt rgb 'red' lw 1.2 pt 7 ps 0.5 title '" << s2.name << "'\n";
  gp.send1d( boost::make_tuple( x1, y1, err1, err2 ) );
  gp.send1d( boost::make_tuple( x1, y1 ) );
  gp.send1d( boost::make_tuple( x2, y2 ) );

  if (b == PDF)
    plot_setup::open_pdf(p);
}

void plot::plot_lines(vector<Line> lines, PlotDetails p, Backend b)
{

  if (lines.size() ==0) return;

  Gnuplot gp;
  plot_setup::setup_gnuplot(gp, b, p);
  gp << "plot ";
  for (int i=0; i<lines.size()-1; ++i) {
    gp << "'-' with lines title '" <<  lines[i].name << "', ";
  }
  gp << "'-' with lines title '" <<  lines.back().name << "'\n";
  for (int i=0; i<lines.size(); ++i) {
    gp.send1d( boost::make_tuple( lines[i].x, lines[i].y ) );
  }

  if (b == PDF)
    plot_setup::open_pdf(p);
}

/*
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
*/
