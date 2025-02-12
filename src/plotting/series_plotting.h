#ifndef SERIES_PLOTTING_H
#define SERIES_PLOTTING_H

#include <vector>
#include <string>
#include <functional>

#include "plot_types.h"

struct Series {
  const std::vector<double>& series;
  std::string name;
};

struct Line {
  std::vector<double> x;
  std::vector<double> y;
  std::string name;
};

struct LineGenerator {
  std::function< double(const std::vector<double>&, unsigned int) > result_gen;
  std::string method_name;
};


namespace plot { 
  void plot_series(Series& s, PlotDetails p);

  void plot_many_series(std::vector<Series>& vs, PlotDetails p);

  void plot_series_diff(Series& s1, Series& s2, PlotDetails p);

  void plot_lines(std::vector<Line> ls, PlotDetails p);

  void plot_lines_generated(const std::vector<double>&, std::vector<unsigned int>, std::vector<LineGenerator>, PlotDetails p);

  void plot_lines_generated_ucr_average(const std::vector<std::string>& dataset_names, std::string dataset_filepath, unsigned int ds_size,
					      std::vector<unsigned int>, std::vector<LineGenerator>, PlotDetails p);

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
}

#endif
