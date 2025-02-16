#ifndef PLOT_TYPES_H
#define PLOT_TYPES_H

#include <gnuplot-iostream.h>
#include <unistd.h>


enum Backend {
  PDF,
  SVG,
  X11,
};

struct PlotDetails {
  std::string title;
  std::string xlabel;
  std::string ylabel;
  std::string filepath;
  Backend b;
};

namespace plot_setup {
  void setup_gnuplot(Gnuplot& gp, PlotDetails p);
  void open_pdf(PlotDetails p);
}


#endif

