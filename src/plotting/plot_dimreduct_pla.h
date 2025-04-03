#ifndef PLOT_DIMREDUCT_PLA_H
#define PLOT_DIMREDUCT_PLA_H

/**
 * @file plot_dimreduct_pla.h
 * @brief Header file containing methods for plotting (Adaptive) Piecewise Linear Approximations
 */

#include "pla.h"
#include <vector>
#include <string>
#include "plot_types.h"

/**
 * @brief plot_pla Is a namespace containing definitions of functions to plot PLA and APLA methods
 */
namespace plot_pla {
  /**
   * @brief plot_pla plots a series with its PLA approximation render on top.
   * @param series Is the time series you want plotted in entirity.
   * @param data_name Is the name of the data passed as series.
   * @param num_params Is the number of parameters allocated to PLA, this is double to the number of segments PLA will use.
   * @param p Is the PlotDetails for the rendered plot.
   */
  void plot_pla(const std::vector<double>& series, std::string data_name, unsigned int num_params, PlotDetails p);
  /**
   * @brief plot_any_apla plots a series over time with its approximation by an adaptive PLA implementation on top.
   * @param series Is the time series you want plotted in entirity.
   * @param data_name Is the name of the data passed as series.
   * @param apla_func Is a function that takes the series you pass in and returns an adaptive PLA representation of it. The form recognised is a tuple of two doubles and unsigned int where the doubles represent a linear function upon the region defined by previous index to the current index (the unsigned integer).
   * @param drt_name Is the name of the dimension reduction technique passed eg. SAPLA.
   * @param p Is the PlotDetails for the rendered plot.
   */
  void plot_any_apla(const std::vector<double>& series, std::string data_name
			    , std::function<const Seqddt(const Seqd&)> apla_func
			    , std::string drt_name
			    , PlotDetails p);
}

#endif
