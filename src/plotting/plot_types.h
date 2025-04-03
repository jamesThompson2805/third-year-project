#ifndef PLOT_TYPES_H
#define PLOT_TYPES_H

/**
 * @file plot_types.h
 * @brief Header file for simplifying actions using gnuplot-iostream.
 */

#include <gnuplot-iostream.h>


/**
 * @brief Enum containing various rendering options for graph.
 * Options are:
 * PDF to render as pdf at <filepath>/<title>.pdf
 * SVG to render as SVG at <filepath>/title.svg
 * X11 to render as X11 window
 */
enum Backend {
  PDF,
  SVG,
  X11,
};

/**
 * @brief Struct containing parameters to be passed to any function wanting to render a plot.
 * @param title is the title of the plot, appearing at top of graph and in the filename.
 * @param xlabel is the label of the x-axis.
 * @param ylabel is the label of the y-axis.
 * @param filepath is the relative or absolute filepath for the image to be saved at, for X11 this is ignored.
 * @param b is the Backend enum chosen to render as.
 */
struct PlotDetails {
  std::string title;
  std::string xlabel;
  std::string ylabel;
  std::string filepath;
  Backend b;
};

/**
 * @brief plot_setup is a namespace containing utilities for setting up gnuplot to render the plot.
 */
namespace plot_setup {
  /**
   * @brief setup_gnuplot should be passed an initialised instance of Gnuplot to configure ready to create a plot.
   * The function uses the details in the PlotDetails to configure the metadata such as size, output, font, filepath and axis labels.
   * @param gp Reference to instance of Gnuplot from gnuplot-iostream libary.
   * @param p PlotDetails instance passed in to configure gnuplot for graph plotting.
   */
  void setup_gnuplot(Gnuplot& gp, PlotDetails p);
  /**
   * @brief open_pdf should be called after graph creation and with same instance of PlotDetails to open specifically pdf plot using firefox.
   * @param p Is PlotDetails instance that should be the same as passed to setup_gnuplot that should have been called before addition of graph.
   */
  void open_pdf(PlotDetails p);
}


#endif

