/**
 * @file demo.cpp
 *
 * @brief This file offers a demo of some of the capabilities of the project designed for the presentation.
 */
#include "random_walk.h"
#include "ucr_parsing.h"
#include "z_norm.h"

#include "bottom_up.h"
#include "pla.h"
#include "double_window.h"
#include "dac_curve_fitting.h"
#include "exact_dp.h"
#include "apla_segment_and_merge.h"
#include "swing.h"
#include "conv_double_window.h"
#include "sliding_window.h"

#include "r_tree.h"
#include "lower_bounds_apla.h"

#include "plotting/series_plotting.h"
#include "plotting/plot_dimreduct_pla.h"
#include "evaluations/capla.h"

#include <iostream>

#include <vector>
#include <string>
using std::vector, std::string;

/**
 * @brief The demo function automatically runs all sections of the demo in order, displaying a full end to end demonstration. It requires floorp browser to run.
 */
void demo()
{
  std::cout << "Third Year Project Demonstration" << std::endl;

  std::cout << " > hit enter to continue" << std::endl;
  string line;
  std::getline(std::cin, line);

  /**** LOAD REAL DATA ****************/
  std::cout << "Load UCR Dataset" << std::endl;
  string ucr_datasets_loc = "external/data/UCRArchive_2018/";
  vector<string> datasets = ucr_parsing::parse_folder_names(ucr_datasets_loc);
  std::cout << std::endl;

  unsigned int di = 30;
  std::cout << " > Dataset ID :  ";
  std::cin >> di;

  vector<double> real_data = parse_ucr_dataset(datasets[di], ucr_datasets_loc,  ucr_parsing::DatasetType::TRAIN_APPEND_TEST);
  vector<double> short_real( real_data.begin(), real_data.begin()+2000);
  z_norm::z_normalise(short_real);

  std::cout << " > Size of loaded dataset : " << real_data.size() << std::endl; 
  Series s_real = { real_data, datasets[di] };
  Series s_short_real = { short_real, datasets[di] };
  PlotDetails p_real = { "Plot of " + datasets[di], "Time", "Value", "./img/demo/", PDF };
  PlotDetails p_short_real = { "Plot of subsequence of " + datasets[di], "Time", "Value", "./img/demo/", PDF };
  //plot::plot_series(s_real, p_real);
  plot::plot_series(s_short_real, p_short_real);

  std::cout << " > hit enter to continue" << std::endl;
  std::getline(std::cin, line);
  std::getline(std::cin, line);

  /***** TEST SOME DRT'S ************/
  std::cout << "Display some Dimension Reduction Techniques on dataset" << std::endl;

  auto d_w_apla_f = [](const vector<double>& s, unsigned int num_params){ return pla::apla_to_seq(d_w::simple_pla(s, num_params, 5, 5)); };
  auto d_w_proj_apla_f = [](const vector<double>& s, unsigned int num_params){ return pla::apla_to_seq(d_w::y_proj_pla(s, num_params, 5, 5)); };
  auto d_w_proj_apla_f_uncompr = [](const vector<double>& s, unsigned int num_params){ return d_w::y_proj_pla(s, num_params, 5, 5); };

  auto exact_apla_f = [](const vector<double>& s, unsigned int num_params){ return pla::apla_to_seq(exact_dp::min_l2_pla(s, num_params)); }; 

  auto rdp_f_uncompr = [&](const Seqd& s, unsigned int parameter){ 
    auto apla = dac_curve_fitting::dac_linear(s, 0.1);
    if (apla.size() < parameter/3) segmerge::segment_to_dim(s,apla,parameter);
    if (apla.size() > parameter/3) segmerge::merge_to_dim(s,apla,parameter);
    return apla;
  };
  auto rdp_f = [&](const Seqd& s, unsigned int parameter) { return pla::apla_to_seq(rdp_f_uncompr(s,parameter)); };
  auto bottom_up_f_uncompr = [&](const Seqd& s, unsigned int parameter){ 
    auto apla = bottom_up::bottom_up(s, 0.1, bottom_up::se);
    if (apla.size() < parameter/3) segmerge::segment_to_dim(s,apla,parameter);
    if (apla.size() > parameter/3) segmerge::merge_to_dim(s,apla,parameter);
    return apla;
  };
  auto bottom_up_f = [&](const Seqd& s, unsigned int parameter) { return pla::apla_to_seq(bottom_up_f_uncompr(s,parameter)); };
  auto sw_f_uncompr = [&](const Seqd& s, unsigned int parameter){ 
    auto apla = sw::sliding_window(s, 0.1);
    if (apla.size() < parameter/3) segmerge::segment_to_dim(s,apla,parameter);
    if (apla.size() > parameter/3) segmerge::merge_to_dim(s,apla,parameter);
    return apla;
  };
  auto sw_f_compr = [&](const Seqd& s, unsigned int parameter) { return pla::apla_to_seq(sw_f_uncompr(s,parameter)); };
  auto swing_f_uncompr = [&](const Seqd& s, unsigned int parameter){ 
    auto apla = swing::swing(s, 0.1);
    if (apla.size() < parameter/3) segmerge::segment_to_dim(s,apla,parameter);
    if (apla.size() > parameter/3) segmerge::merge_to_dim(s,apla,parameter);
    return apla;
  };
  auto swing_f_compr = [&](const Seqd& s, unsigned int parameter) { return pla::apla_to_seq(swing_f_uncompr(s,parameter)); };

  auto sw_f = [&](const Seqd& s, unsigned int parameter) { return pla::apla_to_seq(sw_f_uncompr(s,parameter)); };
  auto swing_f = [&](const Seqd& s) { return swing::swing(s,0.1); };



  unsigned int input;
  std::cout << " 1 : Average Value Double Window	 2 : Interval Projection Double Window" << std::endl;
  std::cout << " 3 : Exact Dynamic Programming		 4 : Top Down" << std::endl;
  std::cout << " 5 : Bottom Up				 6 : Sliding Window" << std::endl;
  std::cout << " 7 : SWING Filter			 8 : Show Evaluations" << std::endl;
  std::cout << " 0 : Exit Loop" << std::endl;
  std::cout << " > Choose an adaptive PLA method: "; 
  std::cin >> input;
  while (input != 0) {
    if (input < 1 || input > 8)
      break;

    if (input == 8) {
      auto d_w_apla_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s, parameter, d_w_apla_f); };
      LineGenerator d_w_apla_gen_l2 = { d_w_apla_gen_l2_f, "AV PLA" };
      auto d_w_proj_apla_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s, parameter, d_w_proj_apla_f); };
      LineGenerator d_w_proj_apla_gen_l2 = { d_w_proj_apla_gen_l2_f, "IP PLA" };
      auto rdp_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s,parameter,rdp_f); };
      LineGenerator rdp_gen_l2 = { rdp_gen_l2_f, "RDP" };
      auto bot_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s,parameter,bottom_up_f); };
      LineGenerator bot_gen_l2 = { bot_gen_l2_f, "B-U" };
      auto sw_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s,parameter,sw_f_compr); };
      LineGenerator sw_gen_l2 = { sw_gen_l2_f, "SW" };
      auto swing_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s,parameter,swing_f_compr); };
      LineGenerator swing_gen_l2 = { swing_gen_l2_f, "SWING" };
      vector<LineGenerator> l2_generators = {
	    d_w_apla_gen_l2
	    , d_w_proj_apla_gen_l2
	    //, apla_gen_l2
	    , rdp_gen_l2
	    , bot_gen_l2
	    , sw_gen_l2
	    , swing_gen_l2
      };
      PlotDetails p_l2 = { "Mean Euclidean Distance of DRTs on " + datasets[di], "Target Dimension m", "L2", "img/demo/", PDF };
      std::vector<Seqd> dataset_samples;
      real_data.resize(10'000);
      z_norm::z_normalise(real_data);
      for (int i=0; i<10; i++) {
	vector<double> dataset_sample(real_data.begin()+1000*i, real_data.begin()+1000*i+1000 );
	z_norm::z_normalise(dataset_sample);
	dataset_samples.push_back(dataset_sample);
      }
      plot::plot_mean_bars_generated(dataset_samples, { 45, 90, 180, 270 }, l2_generators, p_l2);
      std::cout << " 1 : Average Value Double Window	 2 : Interval Projection Double Window" << std::endl;
      std::cout << " 3 : Exact Dynamic Programming		 4 : Top Down" << std::endl;
      std::cout << " 5 : Bottom Up				 6 : Sliding Window" << std::endl;
      std::cout << " 7 : SWING Filter			 8 : Show Evaluations" << std::endl;
      std::cout << " 0 : Exit Loop" << std::endl;
      std::cout << " > Choose an adaptive PLA method: "; 
      std::cin.sync();
      std::cin >> input;

      continue;
    }

    unsigned int num_params=100;
    std::cout << " > Choose data storage of approximation: ";
    std::cin >> num_params;

    Seqd approx;
    std::string name;
    switch (input) {
      case 1: 
	approx = d_w_apla_f(short_real, num_params);
	name = "Average Value Double Window";
	break;
      case 2: 
	approx = d_w_proj_apla_f(short_real, num_params);
	name = "Interval Projection Double Window";
	break;
      case 3: 
	approx = exact_apla_f(short_real, num_params);
	name = "Exact Dynamic Programming";
	break;
      case 4: 
	approx = rdp_f(short_real, num_params);
	name = "Top Down";
	break;
      case 5: 
	approx = bottom_up_f(short_real, num_params);
	name = "Bottom Up";
	break;
      case 6: 
	approx = sw_f(short_real, num_params);
	name = "Sliding Window";
	break;
      default:
	approx = swing_f_compr(short_real, num_params);
	name = "SWING Filter";
	break;
    }
    Series s_approx = { approx, name };
    PlotDetails pd_drt = { "Plot of DRT upon "+datasets[di], "Time", "Value", "./img/demo/", PDF };
    std::vector<Series> v2 = { s_short_real, s_approx };
    plot::plot_many_series(v2, pd_drt);

    std::cout << " 1 : Average Value Double Window	 2 : Interval Projection Double Window" << std::endl;
    std::cout << " 3 : Exact Dynamic Programming		 4 : Top Down" << std::endl;
    std::cout << " 5 : Bottom Up				 6 : Sliding Window" << std::endl;
    std::cout << " 7 : SWING Filter			 8 : Show Evaluations" << std::endl;
    std::cout << " 0 : Exit Loop" << std::endl;
    std::cout << " > Choose an adaptive PLA method: "; 
    std::cin.sync();
    std::cin >> input;
  }

  /****** FIND SOME HEART-BEATS *****/
  std::cout << "\nUse similarity search to find heartbeats" << std::endl;
  std::cout << "hit enter to continue: " << std::endl;
  std::getline(std::cin, line);
  std::getline(std::cin, line);

  di = 30;
  real_data = parse_ucr_dataset(datasets[di], ucr_datasets_loc,  ucr_parsing::DatasetType::TEST);
  z_norm::z_normalise(real_data);
  short_real = std::vector<double>( real_data.begin(), real_data.begin()+2000);
  z_norm::z_normalise(short_real);

  const unsigned int seq_size = 50;
  const unsigned int NS = 5;
  auto retrieval_f = [](const unsigned int& i, const vector<double>& q) {
    return std::vector<std::array<const double*,2>>( {{ q.data()+i, q.data()+i+seq_size-1 }} );
  };
  RTree<apla_bounds::AplaMBR<NS>, unsigned int> r_tree(40,10 , apla_bounds::mbr_area<NS> , apla_bounds::mbr_merge<NS> , apla_bounds::dist_to_mbr_sqr<NS>);
  auto vec_of_mbrs = apla_bounds::vec_to_subseq_mbrs<NS>(short_real,seq_size,bottom_up_f_uncompr);
  for (int i=0; i<vec_of_mbrs.size(); i++) {
    r_tree.insert(vec_of_mbrs[i], i);
  }
  Seqd query(short_real.begin() + 51, short_real.begin() + 101);
  Seqd query_x; for (int i=0; i<query.size(); i++) query_x.push_back( 51 + i );
  Seqd short_real_x; for (int i=0; i<short_real.size(); i++) short_real_x.push_back(i);
  std::vector<Line> vl;
  Line short_real_line = { short_real_x, short_real, "ECG Data" };
  vl.push_back(short_real_line);
  auto sim_results = r_tree.sim_search_exact(query, 7.0, retrieval_f, short_real);
  std::cout << "Start Indexes for similar subsequences : " << std::endl;
  for (const auto& [l_ptr,r_ptr] : sim_results) {
    std::cout << l_ptr - short_real.data() << std::endl;
    Seqd res(l_ptr, r_ptr+1);
    Seqd res_x; for (int i=0; i<res.size(); i++) res_x.push_back( l_ptr - short_real.data() + i );
    vl.push_back({res_x, res, "distance 2 away"});
  }
  vl.push_back({query_x, query, "query"});
  std::cout << std::endl;
  PlotDetails p_sim = { "Plot of subsequences close to first heart beat in ECG (epsilon=7.0)", "Time", "Voltage", "./img/demo/", PDF };
  plot::plot_lines(vl,p_sim);
}
