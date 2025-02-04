#ifndef SERIES_PLOTTING_H
#define SERIES_PLOTTING_H

#include <vector>
#include <string>
#include <functional>

struct Series {
  const std::vector<double>& series;
  std::string name;
};

struct Line {
  std::vector<double>& x;
  std::vector<double>& y;
  std::string name;
};

struct PlotLabels {
  std::string title;
  std::string xlabel;
  std::string ylabel;
};

void plot_series(Series& s, std::string filepath);

void plot_many_series(std::vector<Series>& vs);

void plot_series_diff(Series& s1, Series& s2, std::string filepath);

void plot_lines(std::vector<Line> ls, PlotLabels pl);

void plot_ucr_drt_mse(const std::vector<std::string>& datasets
		      , std::string datasets_loc
		      , unsigned int ds_start
		      , unsigned int ds_end
		      , unsigned int max_num_params
		      , std::function<std::vector<double> (const std::vector<double>&, unsigned int num_params)> drt);

void plot_ucr_drt_maxdev(const std::vector<std::string>& datasets
		      , std::string datasets_loc
		      , unsigned int ds_start
		      , unsigned int ds_end
		      , unsigned int max_num_params
		      , std::function<std::vector<double> (const std::vector<double>&, unsigned int num_params)> drt);

#endif
