#include "series_plotting.h"

#include "gnuplot-iostream.h"
#include "plotting_constants.h"
#include <boost/tuple/tuple.hpp>

using std::vector;

void plot_series(Series& s1)
{
  vector<double> x,y;
  int i=0;
  for (const auto& f : s1.series) {
    x.emplace_back( (double) i);
    y.emplace_back( (double) f);
    ++i;
  }
  Gnuplot gp;
  gp << "set terminal x11\n";
  gp << "set term x11 size 1280 800 \n";
  gp << "set xlabel 'time'\n";
  gp << "set ylabel 'series val'\n";
  gp << "set title 'Plot of series " << s1.name << "'\n";
  gp << "plot '-' with lines title '" <<  s1.name << "'\n";
  gp.send1d( boost::make_tuple( x, y ) );
}

void plot_many_series(vector<Series>& vs)
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
  using namespace gp_constants;
  gp << "set terminal " << GNUPLOT_TERMINAL << "\n";
  gp << "set term " << GNUPLOT_TERMINAL <<" size " << GNUPLOT_SIZE_X << " " << GNUPLOT_SIZE_Y << "\n";
  gp << "set xlabel 'time'\n";
  gp << "set ylabel 'series val'\n";
  gp << "set title 'Plot of multiple time series \n";
  gp << "plot ";
  for (int i=0; i<vs.size()-1; ++i) {
    gp << "'-' with lines title '" <<  vs[i].name << "', ";
  }
  gp << "'-' with lines title '" <<  vs.back().name << "'\n";
  for (int i=0; i<vs.size(); ++i) {
    gp.send1d( boost::make_tuple( vx[i], vy[i] ) );
  }
}

int min(int x, int y) { if (x<y) return x; return y;}

void plot_series_diff(Series& s1, Series& s2)
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
  using namespace gp_constants;
  gp << "set terminal " << GNUPLOT_TERMINAL << "\n";
  gp << "set term " << GNUPLOT_TERMINAL <<" size " << GNUPLOT_SIZE_X << " " << GNUPLOT_SIZE_Y << "\n";
  gp << "set xlabel 'time'\n";
  gp << "set ylabel 'series val'\n";
  gp << "set title 'Comparison of series " << s1.name << " and " << s2.name << "'\n";
  gp << "plot '-' dt 3 with yerrorlines title 'difference'"
     << ", '-' with lines title '" << s1.name << "'"
     << ", '-' with lines title '" << s2.name << "'\n";
  gp.send1d( boost::make_tuple( x1, y1, err1, err2 ) );
  gp.send1d( boost::make_tuple( x1, y1 ) );
  gp.send1d( boost::make_tuple( x2, y2 ) );
}

void plot_lines(vector<Line> lines, PlotLabels pl)
{

  if (lines.size() ==0) return;

  Gnuplot gp;
  using namespace gp_constants;
  gp << "set terminal " << GNUPLOT_TERMINAL << "\n";
  gp << "set term " << GNUPLOT_TERMINAL <<" size " << GNUPLOT_SIZE_X << " " << GNUPLOT_SIZE_Y << "\n";
  gp << "set xlabel '" << pl.xlabel << "'\n";
  gp << "set ylabel '" << pl.ylabel << "'\n";
  gp << "set title '" << pl.title << "'\n";
  gp << "plot ";
  for (int i=0; i<lines.size()-1; ++i) {
    gp << "'-' with lines title '" <<  lines[i].name << "', ";
  }
  gp << "'-' with lines title '" <<  lines.back().name << "'\n";
  for (int i=0; i<lines.size(); ++i) {
    gp.send1d( boost::make_tuple( lines[i].x, lines[i].y ) );
  }
}
