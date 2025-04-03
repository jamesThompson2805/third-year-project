#ifndef PLOT_DIMREDUCT_PAA_H
#define PLOT_DIMREDUCT_PAA_H

/**
 * @file plot_dimreduct_paa.h
 * @brief Header file containing methods for plotting (Adaptive) Piecewise Aggregate Approximations
 */

#include <string>
#include <tuple>
#include <functional>

#include "pla.h"
#include "plot_types.h"

/**
 * @brief plot_paa Is a namespace containing definitions of functions to plot PAA and APCA methods
 */
namespace plot_paa {
  /**
   * @brief plot_paa plots a series with its PAA approximation render on top.
   * @param series Is the time series you want plotted in entirity.
   * @param data_name Is the name of the data passed as series.
   * @param num_params Is the number of parameters allocated to PAA, this is identical to the number of segments PAA will use.
   * @param p Is the PlotDetails for the rendered plot.
   */
  void plot_paa(const Seqd& series, std::string data_name, unsigned int num_params, PlotDetails p);
  /**
   * @brief plot_any_apaa plots a series over time with its approximation by an adaptive PAA implementation on top.
   * @param series Is the time series you want plotted in entirity.
   * @param data_name Is the name of the data passed as series.
   * @param to_apaa Is a function that takes the series you pass in and returns an adaptive PAA representation of it. The form recognised is a tuple of double and unsigned int where the double is the value upon the region defined by previous index to the current index (the unsigned integer).
   * @param drt_name Is the name of the dimension reduction technique passed eg. APCA.
   * @param p Is the PlotDetails for the rendered plot.
   */
  void plot_any_apaa(const Seqd& series , std::string data_name , std::function< std::vector<std::tuple<double, unsigned int>>(const Seqd&) > to_apaa , std::string drt_name , PlotDetails p);
};

#endif
