#include "plot_types.h"
/**
 * @file plot_types.cpp
 * @brief Implementations of functions defined in plot_types.h .
 */

/**
 * @brief setup_gnuplot should be passed an initialised instance of Gnuplot to configure ready to create a plot.
 * The function uses the details in the PlotDetails to configure the metadata such as size, output, font, filepath and axis labels.
 * @param gp Reference to instance of Gnuplot from gnuplot-iostream libary.
 * @param p PlotDetails instance passed in to configure gnuplot for graph plotting.
 */
void plot_setup::setup_gnuplot(Gnuplot &gp, PlotDetails p)
{
  double width = 32;
  double height = 19;
  switch (p.b) {
    case PDF:
      gp << "set term pdfcairo enhanced color dashed font 'Verdana, 24' rounded size "<<width<<"cm, "<<height<<"cm\n";
      gp << "set output '" << p.filepath + p.title <<".pdf'\n";
      break;
    case SVG: 
      gp << "set term svg dashed font 'Verdana, 14' rounded\n";
      gp << "set output '" << p.filepath + p.title <<".svg'\n";
      break;
    case X11:
      gp << "set term x11 dashed font 'Verdana, 14' rounded size 32cm, 19.2cm\n";
  }
  gp << "set loadpath './external/libs'\n";
  gp << "load 'Set1Palette.plt'\n";
  gp << "set xlabel '" << p.xlabel <<"'\n";
  gp << "set ylabel '" << p.ylabel << "'\n";
  gp << "set title '" << p.title << "'\n";
}

/**
 * @brief open_pdf should be called after graph creation and with same instance of PlotDetails to open specifically pdf plot using firefox.
 * @param p Is PlotDetails instance that should be the same as passed to setup_gnuplot that should have been called before addition of graph.
 */
void plot_setup::open_pdf(PlotDetails p)
{
  if (p.b == PDF) {
    std::string command = "firefox '" + p.filepath + p.title + ".pdf'";
    int a = system(command.c_str()); 
  }
}
