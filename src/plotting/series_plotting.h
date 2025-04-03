#ifndef SERIES_PLOTTING_H
#define SERIES_PLOTTING_H

#include <vector>
#include <string>
#include <functional>

#include "plot_types.h"

/**
 * @file series_plotting.h
 * @brief Header file containing many functions for plotting series.
 */

/**
 * @brief Struct stores a series and its name together.
 * @param series Provides constant reference to data to be plotted.
 * @param name Is name you want associated with series.
 */
struct Series {
  const std::vector<double>& series;
  std::string name;
};

/**
 * @brief Struct slightly more complicated than Series also specifying the x values for the corresponding y.
 * Functions using this may or may not check that x and y are same length, this should be maintained for usage.
 * @param x Is the x-values of the data.
 * @param y Is the y-values of the data
 * @param name Is the name you want associated with these coordinates.
 */
struct Line {
  std::vector<double> x;
  std::vector<double> y;
  std::string name;
};

/**
 * @brief Struct stores a method for plotting approximations that vary according to some unsigned int parameter.
 * @param result_gen Is a function that provides better or worse approximations of the passed in series based on parameter unsigned int.
 * @param method_name Is the name you want associated with the function.
 */
struct LineGenerator {
  std::function< double(const std::vector<double>&, unsigned int) > result_gen;
  std::string method_name;
};


/**
 * @brief plot namespace defines useful functions for plotting series and actions performed on series.
 */
namespace plot { 
  /**
   * @brief plot_series plots a series.
   * @param s Is the Series to be plotted.
   * @param p Is the PlotDetails.
   */
  void plot_series(Series& s, PlotDetails p);

  /**
   * @brief plot_many_series plots multiple series.
   * @param vs Is the vector of all Series to be plotted.
   * @param p Is the PlotDetails.
   */
  void plot_many_series(std::vector<Series>& vs, PlotDetails p);
  /**
   * @brief barplot_many_series plots multiple series, here each series is viewed as a value for each set of bars.
   * @param vs Is the vector of all Series to be plotted.
   * @param p Is the PlotDetails.
   */
  void barplot_many_series(std::vector<Series>& vs, PlotDetails p);

  /**
   * @brief plot_series_diff plots two series and adds error lines imbetween to highlight their disagreements.
   * @param s1 Is the first Series to be plotted.
   * @param s2 Is the second Series to be plotted.
   * @param p Is the PlotDetails.
   */
  void plot_series_diff(Series& s1, Series& s2, PlotDetails p);

  /**
   * @brief plot_lines plots all of the lines passed to it using the PlotDetails also passed.
   * @param s2 ls the vector of lines to be plotted.
   * @param p Is the PlotDetails.
   * The function will automatically scale to ensure that all lines are accommodated onto the plot.
   */
  void plot_lines(std::vector<Line> ls, PlotDetails p);

  /**
   * @brief plot_lines_generated plots all of the lines passed to it using the PlotDetails also passed.
   * @param s ls the baseline series.
   * @param x Is the underlying x-axis to plot on (and the values passed to the parameter of LineGenerator).
   * @param y_gens Is the vector of LineGenerator that provide method to vary approximation of s.
   * @param p Is the PlotDetails to render the graph.
   */
  void plot_lines_generated(const std::vector<double>& s, std::vector<unsigned int> x, std::vector<LineGenerator> y_gens, PlotDetails p);

  /**
   * @brief plot_lines_generated_ucr_average plots all of the lines passed to it using the PlotDetails also passed, utilising multiple UCR datasets.
   * @param dataset_names Is the vector of names of the UCR datasets used.
   * @param dataset_filepath Is the filepath to the set of the UCR datasets, these can be downloaded from https://www.cs.ucr.edu/%7Eeamonn/time_series_data_2018/.
   * @param ds_size Is the size the datasets will be set to. This should be small enough such that none are below this size beforehand
   * @param x Is the underlying x-axis to plot on (and the values passed to the parameter of LineGenerator).
   * @param y_gens Is the vector of LineGenerator that provide method to vary approximation of s.
   * @param p Is the PlotDetails to render the graph.
   * Multiple datasets are used via ensuring all are z-normalised and hence being capable to take the mean of all of their results, providing the single line.
   */
  void plot_lines_generated_ucr_average(const std::vector<std::string>& dataset_names, std::string dataset_filepath, unsigned int ds_size,
					      std::vector<unsigned int> x, std::vector<LineGenerator> y_gens, PlotDetails p);
}

#endif
