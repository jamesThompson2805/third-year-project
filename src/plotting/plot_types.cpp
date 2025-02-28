#include "plot_types.h"

void plot_setup::setup_gnuplot(Gnuplot &gp, PlotDetails p)
{
  double width = 32;
  double height = 32;
  switch (p.b) {
    case PDF:
      gp << "set term pdfcairo enhanced color dashed font 'Verdana, 30' rounded size "<<width<<"cm, "<<height<<"cm\n";
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

void plot_setup::open_pdf(PlotDetails p)
{
  if (p.b == PDF) {
    std::string command = "firefox '" + p.filepath + p.title + ".pdf'";
    int a = system(command.c_str()); 
  }
}
