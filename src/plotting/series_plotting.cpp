#include "series_plotting.h"

/**
 * @file series_plotting.cpp
 * @brief Implementations of functions defined in series_plotting.h .
 */

#include "gnuplot-iostream.h"
#include "error_measures.h"
#include "plot_types.h"
#include <algorithm>
#include <boost/tuple/tuple.hpp>

#include "ucr_parsing.h"
#include "z_norm.h"

#include "pgbar.hpp"

using std::vector;

void plot::plot_series(Series& s1, PlotDetails p)
{
  vector<double> x,y;
  int i=0;
  for (const auto& f : s1.series) {
    x.emplace_back( (double) i);
    y.emplace_back( (double) f);
    ++i;
  }
  Gnuplot gp;
  plot_setup::setup_gnuplot(gp, p);
  gp << "plot '-' with lines lt rgb 'blue' lw 1.2 title '" <<  s1.name << "'\n";
  gp.send1d( boost::make_tuple( x, y ) );

  plot_setup::open_pdf(p);
}

void plot::plot_many_series(vector<Series>& vs, PlotDetails p)
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
  plot_setup::setup_gnuplot(gp, p);
  gp << "plot ";
  for (int i=0; i<vs.size()-1; ++i) {
    gp << "'-' with lines dt " << i+1 << " title '" <<  vs[i].name << "', ";
  }
  gp << "'-' with lines dt " << vs.size() << " title '" <<  vs.back().name << "'\n";
  for (int i=0; i<vs.size(); ++i) {
    gp.send1d( boost::make_tuple( vx[i], vy[i] ) );
  }

  plot_setup::open_pdf(p);
}

void plot::barplot_many_series(const std::vector<Series> &vs, const std::vector<std::string>& group_labels, PlotDetails p)
{
  if (vs.size() == 0) return;

  Gnuplot gp;
  plot_setup::setup_gnuplot(gp, p);

  gp << "set style line 2 lc rgb 'black' lt 1 lw 1\n";
  gp << "set style data histogram\n";
  gp << "set style histogram cluster gap 1\n";
  gp << "set style fill pattern border -1\n";
  gp << "set boxwidth 0.9\n";
  gp << "set xtics format ''\n";

  //gp << "set key rmargin\n";
  gp << "set key off\n";

  //gp << "set yrange [0:500]\n";

  gp << "set xtics (";
  for (int i=0; i<group_labels.size()-1; i++)
    gp << "'" << group_labels[i] << "' "<<std::to_string(i) <<", ";
  gp << "'" << group_labels.back()
    << "' " << std::to_string(group_labels.size()-1) << ")\n";

  gp << "set grid ytics\n";

  gp << "plot ";
  for (int i=0; i<vs.size()-1; ++i) {
    gp << "'-' title '" <<  vs[i].name << "' ls 2, ";
  }
  gp << "'-' title '" <<  vs.back().name << "' ls 2\n";
  for (int i=0; i<vs.size(); ++i) {
    gp.send1d( boost::make_tuple( vs[i].series ) );
  }
  plot_setup::open_pdf(p);
}

int min(int x, int y) { if (x<y) return x; return y;}

void plot::plot_series_diff(Series& s1, Series& s2, PlotDetails p)
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
  plot_setup::setup_gnuplot(gp, p);
  gp << "plot '-' dt 3 lw 2 with yerrorlines title 'difference'"
     << ", '-' with linespoints lt rgb 'blue' lw 1.2 pt 5 ps 0.5 title '" << s1.name << "'"
     << ", '-' with linespoints lt rgb 'red' lw 1.2 pt 7 ps 0.5 title '" << s2.name << "'\n";
  gp.send1d( boost::make_tuple( x1, y1, err1, err2 ) );
  gp.send1d( boost::make_tuple( x1, y1 ) );
  gp.send1d( boost::make_tuple( x2, y2 ) );

  plot_setup::open_pdf(p);
}

void plot::plot_lines(vector<Line> lines, PlotDetails p)
{

  if (lines.size() ==0) return;

  Gnuplot gp;
  plot_setup::setup_gnuplot(gp, p);
  //gp << "set key off\n";
  gp << "set key outside\n";
  gp << "plot ";
  gp << "'-' with lines lw 4.0 title '" <<  lines[0].name << "', ";
  for (int i=1; i<lines.size()-1; ++i) {
    gp << "'-' with lines lw 4.0 dt " << i+1 << " title '" <<  lines[i].name << "', ";
    //gp << "'-' with lines lw 4.0 title '" <<  lines[i].name << "', ";
  }
  gp << "'-' with lines lw 4.0 dt " << lines.size() << " title '" <<  lines.back().name << "'\n";
  //gp << "'-' with lines lw 4.0 title '" <<  lines.back().name << "'\n";
  for (int i=0; i<lines.size(); ++i) {
    gp.send1d( boost::make_tuple( lines[i].x, lines[i].y ) );
  }

  plot_setup::open_pdf(p);
}

void plot::plot_lines_generated(const std::vector<double>& s, std::vector<unsigned int> x, std::vector<LineGenerator> y_gens, PlotDetails p)
{
  if (y_gens.size() == 0 || x.size() == 0) return;
  vector<double> x_d(x.size());
  std::transform(x.begin(), x.end(), x_d.begin(), [](unsigned int& i){ return (double) i;}); 

  vector<vector<double>> y_vecs;
  vector<Line> lines;
  for (auto& y_f : y_gens) {
    y_vecs.push_back({});
    for (unsigned int i : x) {
      y_vecs.back().push_back( y_f.result_gen(s,i) );
    }
    lines.push_back({x_d, y_vecs.back(), y_f.method_name});
  }
  plot::plot_lines(lines, p);
}
void plot::plot_bars_generated(const std::vector<double>& s, std::vector<unsigned int> x, std::vector<LineGenerator> y_gens, PlotDetails p)
{
  if (y_gens.size() == 0 || x.size() == 0) return;
  vector<std::string> x_d(x.size());
  std::transform(x.begin(), x.end(), x_d.begin(), [](unsigned int& i){ return  std::to_string(i);}); 

  vector<vector<double>> y_vecs;
  vector<Series> lines;
  for (auto& y_f : y_gens) {
    y_vecs.push_back({});
    for (unsigned int i : x) {
      y_vecs.back().push_back( y_f.result_gen(s,i) );
    }
    lines.push_back({ y_vecs.back(), y_f.method_name});
  }
  plot::barplot_many_series(lines, x_d, p);
}
void plot::plot_barsPG_generated(const std::vector<double>& s, std::vector<double> x, std::vector<LinePGGenerator> y_gens, PlotDetails p)
{
  if (y_gens.size() == 0 || x.size() == 0) return;
  vector<std::string> x_d(x.size());
  std::transform(x.begin(), x.end(), x_d.begin(), [](double& i){ return  std::to_string(i);}); 

  vector<vector<double>> y_vecs;
  vector<Series> lines;
  for (auto& y_f : y_gens) {
    y_vecs.push_back({});
    for (double e : x) {
      y_vecs.back().push_back( y_f.result_gen(s,e) );
    }
    lines.push_back({ y_vecs.back(), y_f.method_name});
  }
  plot::barplot_many_series(lines, x_d, p);
}

void plot::plot_lines_generated_ucr_average(const std::vector<std::string>& dataset_names, std::string dataset_filepath, unsigned int ds_size,
				      std::vector<unsigned int> x, std::vector<LineGenerator> y_gens, PlotDetails p)
{
  if (y_gens.size() == 0 || x.size() == 0) return;
  vector<double> x_d(x.size());
  std::transform(x.begin(), x.end(), x_d.begin(), [](unsigned int& i){ return (double) i;}); 

  pgbar::ProgressBar<> bar { pgbar::option::Remains( "-" ),
                             pgbar::option::Filler( "=" ),
                             pgbar::option::Style( pgbar::config::CharBar::Entire ),
                             pgbar::option::RemainsColor( "#A52A2A" ),
                             pgbar::option::FillerColor( 0x0099FF ),
                             pgbar::option::InfoColor( pgbar::color::Yellow ),
                             pgbar::option::Tasks( x.size() * dataset_names.size() ) };
  vector<vector<double>> y_vecs(y_gens.size());
  vector<Line> lines;
  for (unsigned int xi : x) {
    std::for_each(y_vecs.begin(), y_vecs.end(), [](auto& v){ v.push_back(.0); });

    int num_used_datasets=dataset_names.size();
    for (const std::string& d_name : dataset_names) {
      bar.tick();
      // construct the dataset
      vector<double> dataset = ucr_parsing::parse_ucr_dataset(d_name, dataset_filepath, ucr_parsing::DatasetType::TRAIN);
      if (ds_size > dataset.size()) {
	num_used_datasets--;
	continue;
      }

      if (ds_size > 0)
	dataset.resize( std::min( (unsigned int) dataset.size(),ds_size) );
      z_norm::z_normalise(dataset);

      for (int yi=0; yi<y_gens.size(); yi++) {
	auto& y_f = y_gens[yi];
	y_vecs[yi].back() += y_f.result_gen(dataset,xi);
      }
    }
    std::for_each(y_vecs.begin(), y_vecs.end(), [&num_used_datasets](auto& v){ v.back() /= (double) num_used_datasets; });
  }

  for (int yi=0; yi<y_gens.size(); yi++) {
    auto& y_vec = y_vecs[yi];
    lines.push_back( {x_d, y_vec, y_gens[yi].method_name} );
  }
  plot::plot_lines(lines, p);
}

