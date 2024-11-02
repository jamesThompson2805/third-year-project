#ifndef SERIES_PLOTTING_H
#define SERIES_PLOTTING_H

#include <vector>
#include <string>

struct Series {
  const std::vector<float>& series;
  std::string name;
};

struct Line {
  std::vector<float>& x;
  std::vector<float>& y;
  std::string name;
};

struct PlotLabels {
  std::string title;
  std::string xlabel;
  std::string ylabel;
};

void plot_series(Series& s);

void plot_many_series(std::vector<Series>& vs);

void plot_series_diff(Series& s1, Series& s2);

void plot_lines(std::vector<Line> ls, PlotLabels pl);

#endif
