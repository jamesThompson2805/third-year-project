/**
 * @file main.cpp
 *
 * @brief Main implementation highlighting usage of utilities and DRTs
 */
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <tuple>
#include <algorithm>

#include "plotting/plot_dimreduct_pla.h"
#include "ucr_parsing.h"
#include "z_norm.h"
#include "error_measures.h"

#include "paa.h"
#include "apca.h"

#include "pla.h"
#include "double_window.h"
#include "dac_curve_fitting.h"
#include "bottom_up.h"
#include "exact_dp.h"
#include "apla_segment_and_merge.h"
#include "swing.h"
#include "conv_double_window.h"
#include "sliding_window.h"

#include "plotting/series_plotting.h"
#include "plotting/plot_dimreduct_paa.h"

#include "evaluations/general.h"
#include "evaluations/capla.h"
#include "evaluations/e_guarantee_eval.h"

#include "random_walk.h"

#include "r_tree.h"
#include "lower_bounds_apla.h"
#include "sequential_scan.h"

#include <chrono>

#include "demo.cpp"


using std::vector, std::string, std::tuple;

/**
 * @brief The entry point of the project and a testing ground for all the project implements
 *
 * The main function is a large testing box of many of the features implemented in the project, left as usage examples for navigating.
 * It highlights real and synthetic data, applying DRT's including varying parameters and the many plotting capabilities as well.
 */
int main()
{
  using namespace ucr_parsing;

  string ucr_datasets_loc = "external/data/UCRArchive_2018/";
  vector<string> datasets = parse_folder_names(ucr_datasets_loc);

  /*
  for (int i=0; i<datasets.size(); i++) {
    vector<double> dataset = parse_ucr_dataset(datasets[i], ucr_datasets_loc,  DatasetType::TRAIN);
    std::cout << i << " : " << datasets[i] << " : " << dataset.size() << std::endl;
  }
  */

  // Used in reports : 5 is Arrowhead, 13 is Chlorine, 28 is ECG200, 39 is Fifty Words, 128 is yoga, 46 is GuestureMidAirD1
  // 109 is Strawberry, 25 is Dodgers loop day
  unsigned int di = 5;
  vector<double> dataset = parse_ucr_dataset(datasets[di], ucr_datasets_loc,  DatasetType::TRAIN_APPEND_TEST);
  //std::cout << dataset.size() << std::endl;
  z_norm::z_normalise(dataset);

  Series s_ucr = { dataset, "Chlorine Dataset" };
  PlotDetails pd_ucr = { "Plot of ucr rough data", "Time", "Value", "img/walks/", X11 };
  //plot::plot_series(s_ucr, pd_ucr);

  auto paa_f = [](const vector<double>& s, unsigned int num_params){ return paa::paa_to_seq( paa::paa(s, num_params), (s.size() / num_params) + (s.size() % num_params != 0)); };
  auto pla_f = [](const vector<double>& s, unsigned int num_params){ return pla::pla_to_seq(pla::pla(s, num_params), (2*s.size() / num_params) + (2*s.size() % num_params != 0)); };

  auto d_w_apca_f = [](const vector<double>& s, unsigned int num_params){ return d_w::simple_paa(s, num_params, 5, 5); };
  auto d_w_proj_apca_f = [](const vector<double>& s, unsigned int num_params){ return d_w::y_proj_paa(s, num_params, 5, 5); };

  auto d_w_apla_f = [](const vector<double>& s, unsigned int num_params){ return pla::apla_to_seq(d_w::simple_pla(s, num_params, 5, 5)); };
  auto d_w_proj_apla_f = [](const vector<double>& s, unsigned int num_params){ return pla::apla_to_seq(d_w::y_proj_pla(s, num_params, 5, 5)); };
  auto d_w_proj_apla_f_uncompr = [](const vector<double>& s, unsigned int num_params){ return d_w::y_proj_pla(s, num_params, 5, 5); };

  auto exact_apaa_f = [](const vector<double>& s, unsigned int num_params){ return paa::apca_to_seq(exact_dp::min_mse_paa(s, num_params)); }; 
  auto exact_apla_f = [](const vector<double>& s, unsigned int num_params){ return pla::apla_to_seq(exact_dp::min_mse_pla(s, num_params)); }; 

  auto apca_f = [](const vector<double>& s, unsigned int num_params){ return paa::apca_to_seq(apca::apca(s, num_params)); };

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

  vector<double> ldist = { 1.0/3.0, 1.0/3.0, 1.0/3.0};
  vector<double> rdist = { 0.0, 1.0/2.0, 1.0/2.0};
  auto conv_apla_f = [&ldist, &rdist](const vector<double>& s, unsigned int num_params){ return c_d_w::conv_pla(s, num_params, ldist, rdist); }; 

  RandomWalk walk( NormalFunctor(1) ); 
  walk.gen_steps(120);
  //walk.save_walk("./tsv/testwalk1.tsv");
  vector<double> dataset2 = parse_tsv("./tsv/testwalk1.tsv",-1);
  z_norm::z_normalise(dataset2);
  Series s = { dataset2, "Normal Walk" };
  
  /************************************* DEMO ******************************************/
  
  demo();
  return 0;
  
  /*
  dataset.resize(100);
  auto sw_f_2 = [&](const Seqd& s) { return sw_f_uncompr(s,42); };
  plot_pla::plot_any_apla(dataset, "chlorine", sw_f_2, "sliding window", pd_ucr);
  */

  
  /************************** Plot of DRT vs original **************************************
  PlotDetails p = { "RDP nonsense", "Time", "", "img/", PDF };
  auto apla_uncompr = [](const Seqd& s){
    auto apla = dac_curve_fitting::dac_linear(s, 0.1);
    segmerge::merge_to_dim(s, apla, 45);
    return apla;
  };
  plot_any_apla_subseq(dataset, "Arrowhead Scan", apla_uncompr, 0, 999, "RDP", p);
  ************/


  auto capla_mean_f =  capla_eval::generate_mean_DRT(5);
  auto capla_tri_f =  capla_eval::generate_tri_DRT(5);
  auto capla_mean_s1_f =  capla_eval::generate_mean_skip_one_DRT(5);
  auto capla_tri_s1_f =  capla_eval::generate_tri_skip_one_DRT(5);

  /************************** Euclidean Evaluation against parameters *************************/
  // PAA
  auto paa_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s, parameter, paa_f); };
  LineGenerator paa_gen_l2 = { paa_gen_l2_f, "PAA" };
  // APCA
  auto apca_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s, parameter, apca_f); };
  LineGenerator apca_gen_l2 = { apca_gen_l2_f, "APCA" };
  // PLA
  auto pla_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s, parameter, pla_f); };
  LineGenerator pla_gen_l2 = { pla_gen_l2_f, "PLA" };
  // Double Window PLA and Proj
  auto d_w_apla_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s, parameter, d_w_apla_f); };
  LineGenerator d_w_apla_gen_l2 = { d_w_apla_gen_l2_f, "AV PLA" };
  auto d_w_proj_apla_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s, parameter, d_w_proj_apla_f); };
  LineGenerator d_w_proj_apla_gen_l2 = { d_w_proj_apla_gen_l2_f, "IP PLA" };
  // CAPLA hypotheses
  auto c_d_w_mean_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s, parameter, capla_mean_f); };
  LineGenerator c_d_w_mean_gen_l2 = { c_d_w_mean_gen_l2_f, "Faster Double Window PLA(5,5)" };
  auto c_d_w_tri_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s, parameter, capla_tri_f); };
  LineGenerator c_d_w_tri_gen_l2 = { c_d_w_tri_gen_l2_f, "IncrWAV PLA" };
  auto c_d_w_mean_s1_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s, parameter, capla_mean_s1_f); };
  LineGenerator c_d_w_mean_s1_gen_l2 = { c_d_w_mean_s1_gen_l2_f, "SkipWAV PLA" };
  auto c_d_w_tri_s1_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s, parameter, capla_tri_s1_f); };
  LineGenerator c_d_w_tri_s1_gen_l2 = { c_d_w_tri_s1_gen_l2_f, "Tri s1 double window pla" };
  auto apla_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s, parameter, exact_apla_f); };
  LineGenerator apla_gen_l2 = { apla_gen_l2_f, "APLA" };
  auto rdp_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s,parameter,rdp_f); };
  LineGenerator rdp_gen_l2 = { rdp_gen_l2_f, "RDP" };
  auto bot_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s,parameter,bottom_up_f); };
  LineGenerator bot_gen_l2 = { bot_gen_l2_f, "B-U" };

  auto sw_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s,parameter,sw_f_compr); };
  LineGenerator sw_gen_l2 = { sw_gen_l2_f, "SW" };
  auto swing_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s,parameter,swing_f_compr); };
  LineGenerator swing_gen_l2 = { swing_gen_l2_f, "SWING" };
  vector<LineGenerator> l2_generators = {
	paa_gen_l2
	, pla_gen_l2
	, apca_gen_l2
	, d_w_apla_gen_l2
	, d_w_proj_apla_gen_l2
	, c_d_w_mean_s1_gen_l2
	, c_d_w_tri_gen_l2
	//, apla_gen_l2
	, rdp_gen_l2
	, bot_gen_l2
	//, sw_gen_l2
	, swing_gen_l2
  };
  /**************/

  /************************** maxdev Evaluation against parameters ****************************/
  // PAA
  auto paa_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, paa_f); };
  LineGenerator paa_gen_maxdev = { paa_gen_maxdev_f, "PAA" };
  // APCA
  auto apca_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, apca_f); };
  LineGenerator apca_gen_maxdev = { apca_gen_maxdev_f, "APCA" };
  // PLA
  auto pla_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, pla_f); };
  LineGenerator pla_gen_maxdev = { pla_gen_maxdev_f, "PLA" };
  // Double Window PLA and Proj
  auto d_w_apla_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, d_w_apla_f); };
  LineGenerator d_w_apla_gen_maxdev = { d_w_apla_gen_maxdev_f, "AV PLA" };
  auto d_w_proj_apla_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, d_w_proj_apla_f); };
  LineGenerator d_w_proj_apla_gen_maxdev = { d_w_proj_apla_gen_maxdev_f, "IP PLA" };
  // CAPLA hypotheses
  auto c_d_w_mean_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, capla_mean_f); };
  LineGenerator c_d_w_mean_gen_maxdev = { c_d_w_mean_gen_maxdev_f, "Faster Double Window PLA(5,5)" };
  auto c_d_w_tri_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, capla_tri_f); };
  LineGenerator c_d_w_tri_gen_maxdev = { c_d_w_tri_gen_maxdev_f, "IncrWAV PLA" };
  auto c_d_w_mean_s1_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, capla_mean_s1_f); };
  LineGenerator c_d_w_mean_s1_gen_maxdev = { c_d_w_mean_s1_gen_maxdev_f, "SkipWAV PLA" };
  auto c_d_w_tri_s1_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, capla_tri_s1_f); };
  LineGenerator c_d_w_tri_s1_gen_maxdev = { c_d_w_tri_s1_gen_maxdev_f, "Tri s1 double window pla(5,5)" };
  auto apla_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, exact_apla_f); };
  LineGenerator apla_gen_maxdev = { apla_gen_maxdev_f, "APLA" };
  auto rdp_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, rdp_f); };
  LineGenerator rdp_gen_maxdev = { rdp_gen_maxdev_f, "RDP" };
  auto bot_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s,parameter,bottom_up_f); };
  LineGenerator bot_gen_maxdev = { bot_gen_maxdev_f, "B-U" };

  auto sw_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s,parameter,sw_f_compr); };
  LineGenerator sw_gen_maxdev = { sw_gen_maxdev_f, "SW" };
  auto swing_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s,parameter,swing_f_compr); };
  LineGenerator swing_gen_maxdev = { swing_gen_maxdev_f, "SWING" };
  vector<LineGenerator> maxdev_generators = {
	paa_gen_maxdev
	, pla_gen_maxdev
	, apca_gen_maxdev
	, d_w_apla_gen_maxdev
	, d_w_proj_apla_gen_maxdev
	, c_d_w_mean_s1_gen_maxdev
	, c_d_w_tri_gen_maxdev
	//, apla_gen_maxdev
	, rdp_gen_maxdev
	, bot_gen_maxdev
	//, sw_gen_maxdev
	, swing_gen_maxdev
  };
  /**************/

  /************************** Time Evaluation against parameters ******************************/
  // PAA
  auto paa_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, paa_f); };
  LineGenerator paa_gen_time = { paa_gen_time_f, "PAA" };
  // APCA
  auto apca_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, apca_f); };
  LineGenerator apca_gen_time = { apca_gen_time_f, "APCA" };
  // PLA
  auto pla_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, pla_f); };
  LineGenerator pla_gen_time = { pla_gen_time_f, "PLA" };
  // Double Window PLA and Proj
  auto d_w_apla_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, d_w_apla_f); };
  LineGenerator d_w_apla_gen_time = { d_w_apla_gen_time_f, "AV PLA" };
  auto d_w_proj_apla_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, d_w_proj_apla_f); };
  LineGenerator d_w_proj_apla_gen_time = { d_w_proj_apla_gen_time_f, "IP PLA" };
  // CAPLA hypotheses
  auto c_d_w_mean_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, capla_mean_f); };
  LineGenerator c_d_w_mean_gen_time = { c_d_w_mean_gen_time_f, "Faster Double Window PLA" };
  auto c_d_w_tri_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, capla_tri_f); };
  LineGenerator c_d_w_tri_gen_time = { c_d_w_tri_gen_time_f, "IncrWAV PLA" };
  auto c_d_w_mean_s1_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, capla_mean_s1_f); };
  LineGenerator c_d_w_mean_s1_gen_time = { c_d_w_mean_s1_gen_time_f, "SkipWAV PLA" };
  auto c_d_w_tri_s1_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, capla_tri_s1_f); };
  LineGenerator c_d_w_tri_s1_gen_time = { c_d_w_tri_s1_gen_time_f, "Tri s1 double window pla" };
  auto apla_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, exact_apla_f); };
  LineGenerator apla_gen_time = { apla_gen_time_f, "APLA" };
  auto rdp_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, rdp_f); };
  LineGenerator rdp_gen_time = { rdp_gen_time_f, "RDP" };
  auto bot_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s,parameter,bottom_up_f); };
  LineGenerator bot_gen_time = { bot_gen_time_f, "B-U" };

  auto sw_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s,parameter,sw_f_compr); };
  LineGenerator sw_gen_time = { sw_gen_time_f, "SW" };
  auto swing_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s,parameter,swing_f_compr); };
  LineGenerator swing_gen_time = { swing_gen_time_f, "SWING" };
  vector<LineGenerator> time_generators = {
	paa_gen_time
	, apca_gen_time
	, pla_gen_time
	, d_w_apla_gen_time
	, d_w_proj_apla_gen_time
	, c_d_w_mean_s1_gen_time
	, c_d_w_tri_gen_time
	//, apla_gen_time
	, rdp_gen_time
	, bot_gen_time
	//, sw_gen_time
	, swing_gen_time
  };
  /**************/

  /************************** Time Evaluation against size ************************************/
  unsigned int tdim = 150;
  // PAA
  auto paa_gen_size_time_f = [&](const Seqd& s, unsigned int parameter){ 
    Seqd s1(s.begin(),s.begin()+parameter);
    return general_eval::cputime_ms_of_method(s1, tdim, paa_f);
  };
  LineGenerator paa_gen_size_time = { paa_gen_size_time_f, "PAA" };
  // APCA
  auto apca_gen_size_time_f = [&](const Seqd& s, unsigned int parameter){ 
    Seqd s1(s.begin(),s.begin()+parameter);
    return general_eval::cputime_ms_of_method(s1, tdim, apca_f);
  };
  LineGenerator apca_gen_size_time = { apca_gen_size_time_f, "APCA" };
  // PLA
  auto pla_gen_size_time_f = [&](const Seqd& s, unsigned int parameter){ 
    Seqd s1(s.begin(),s.begin()+parameter);
    return general_eval::cputime_ms_of_method(s1, tdim, pla_f);
  };
  LineGenerator pla_gen_size_time = { pla_gen_size_time_f, "PLA" };
  // Double Window PLA and Proj
  auto d_w_apla_gen_size_time_f = [&](const Seqd& s, unsigned int parameter){ 
    Seqd s1(s.begin(),s.begin()+parameter);
    return general_eval::cputime_ms_of_method(s1, tdim, d_w_apla_f);
  };
  LineGenerator d_w_apla_gen_size_time = { d_w_apla_gen_size_time_f, "Double Window PLA" };
  auto d_w_proj_apla_gen_size_time_f = [&](const Seqd& s, unsigned int parameter){ 
    Seqd s1(s.begin(),s.begin()+parameter);
    return general_eval::cputime_ms_of_method(s1, tdim, d_w_proj_apla_f);
  };
  LineGenerator d_w_proj_apla_gen_size_time = { d_w_proj_apla_gen_size_time_f, "Projection PLA" };
  // CAPLA hypotheses
  auto c_d_w_mean_gen_size_time_f = [&](const Seqd& s, unsigned int parameter){ 
    Seqd s1(s.begin(),s.begin()+parameter);
    return general_eval::cputime_ms_of_method(s1, tdim, capla_mean_f);
  };
  LineGenerator c_d_w_mean_gen_size_time = { c_d_w_mean_gen_size_time_f, "Faster Double Window PLA" };
  auto c_d_w_tri_gen_size_time_f = [&](const Seqd& s, unsigned int parameter){ 
    Seqd s1(s.begin(),s.begin()+parameter);
    return general_eval::cputime_ms_of_method(s1, tdim, capla_tri_f);
  };
  LineGenerator c_d_w_tri_gen_size_time = { c_d_w_tri_gen_size_time_f, "Triangular Double Window PLA" };
  auto c_d_w_mean_s1_gen_size_time_f = [&](const Seqd& s, unsigned int parameter){ 
    Seqd s1(s.begin(),s.begin()+parameter);
    return general_eval::cputime_ms_of_method(s1, tdim, capla_mean_s1_f);
  };
  LineGenerator c_d_w_mean_s1_gen_size_time = { c_d_w_mean_s1_gen_size_time_f, "mean s1 double window pla" };
  auto c_d_w_tri_s1_gen_size_time_f = [&](const Seqd& s, unsigned int parameter){ 
    Seqd s1(s.begin(),s.begin()+parameter);
    return general_eval::cputime_ms_of_method(s1, tdim, capla_tri_s1_f);
  };
  LineGenerator c_d_w_tri_s1_gen_size_time = { c_d_w_tri_s1_gen_size_time_f, "Tri s1 double window pla" };
  auto apla_gen_size_time_f = [&](const Seqd& s, unsigned int parameter){ 
    Seqd s1(s.begin(),s.begin()+parameter);
    return general_eval::cputime_ms_of_method(s1, tdim, exact_apla_f);
  };
  LineGenerator apla_gen_size_time = { apla_gen_size_time_f, "DP Exact PLA" };
  auto rdp_gen_size_time_f = [&](const Seqd& s, unsigned int parameter){ 
    Seqd s1(s.begin(),s.begin()+parameter);
    return general_eval::cputime_ms_of_method(s1, tdim, rdp_f);
  };
  LineGenerator rdp_gen_size_time = { rdp_gen_size_time_f, "RDP" };
  auto bot_gen_size_time_f = [&](const Seqd& s, unsigned int parameter){ 
    Seqd s1(s.begin(),s.begin()+parameter);
    return general_eval::cputime_ms_of_method(s1, tdim, bottom_up_f);
  };
  LineGenerator bot_gen_size_time = { bot_gen_size_time_f, "B-U" };
  vector<LineGenerator> time_size_generators = {
	paa_gen_size_time
	, pla_gen_size_time
	, apca_gen_size_time
	, d_w_apla_gen_size_time
	, d_w_proj_apla_gen_size_time
	, c_d_w_mean_gen_size_time
	, c_d_w_tri_gen_size_time
	//, apla_gen_size_time
	//, rdp_gen_size_time
	//, bot_gen_size_time
  };
  /**************/

  /************************** Compression Ratio against epsilon ************************************/
  auto bottom_up = [&](const Seqd& s, double e) { return bottom_up::bottom_up(s, e, bottom_up::maxdev); };

  auto top_down_compr_ratio  = [&](const Seqd& s, double e) { return precision_eval::compr_ratio_of_method(s,e,dac_curve_fitting::dac_linear); };
  auto bottom_up_compr_ratio = [&](const Seqd& s, double e) { return precision_eval::compr_ratio_of_method(s,e,bottom_up); };
  auto sliding_w_compr_ratio = [&](const Seqd& s, double e) { return precision_eval::compr_ratio_of_method(s,e,sw::sliding_window); };
  auto swing_compr_ratio     = [&](const Seqd& s, double e) { return precision_eval::compr_ratio_of_method(s,e,swing::swing); };

  auto top_down_time  = [&](const Seqd& s, double e) { return precision_eval::cputime_of_method(s,e,dac_curve_fitting::dac_linear); };
  auto bottom_up_time = [&](const Seqd& s, double e) { return precision_eval::cputime_of_method(s,e,bottom_up); };
  auto sliding_w_time = [&](const Seqd& s, double e) { return precision_eval::cputime_of_method(s,e,sw::sliding_window); };
  auto swing_time     = [&](const Seqd& s, double e) { return precision_eval::cputime_of_method(s,e,swing::swing); };

  LinePGGenerator top_down_gen = { top_down_compr_ratio, "Top Down" };
  LinePGGenerator bottom_up_gen = { bottom_up_compr_ratio, "Bottom Up" };
  LinePGGenerator sw_gen = { sliding_w_compr_ratio, "Sliding Window" };
  LinePGGenerator swing_gen = { swing_compr_ratio, "Swing" };

  LinePGGenerator top_down_gen_time = { top_down_time, "Top Down" };
  LinePGGenerator bottom_up_gen_time = { bottom_up_time, "Bottom Up" };
  LinePGGenerator sw_eps_gen_time = { sliding_w_time, "Sliding Window" };
  LinePGGenerator swing_eps_gen_time = { swing_time, "Swing" };

  vector<LinePGGenerator> compr_generators = { top_down_gen, bottom_up_gen, sw_gen, swing_gen };
  vector<LinePGGenerator> etime_generators = { top_down_gen_time, bottom_up_gen_time, sw_eps_gen_time, swing_eps_gen_time };

  /*
  vector<unsigned int> compr_ds = { 5, 13, 25, 109 };
  for (unsigned int compr_i=0; compr_i < datasets.size(); compr_i++) {
    dataset = parse_ucr_dataset(datasets[compr_i], ucr_datasets_loc,  DatasetType::TRAIN_APPEND_TEST);
    //std::cout << dataset.size() << std::endl;
    if (dataset.size() < 2000) continue;
    dataset.resize(2000);
    z_norm::z_normalise(dataset);
    std::cout << compr_i <<std::endl;
    if (bottom_up_compr_ratio(dataset, 0.15) < sliding_w_compr_ratio(dataset, 0.15) ) {
      std::cout <<"	" <<compr_i << std::endl;
    }

  }
  */
  

  /**************/




  PlotDetails p_l2 = { "Euclidean Distance of DRTs on Chlorine dataset", "Target Dimension m", "L2", "img/drt_comparisons/", PDF };
  PlotDetails p_maxdev = { "Maximum Deviation of DRTs on Chlorine dataset", "Target Dimension m", "Maximum Deviation", "img/drt_comparisons/", PDF };
  PlotDetails p_time = { "CPU Time of DRTs on Synthetic dataset", "Target Dimension m", "CPU Execution Time (ns)", "img/drt_comparisons/", PDF };
  PlotDetails p_time_size = { "CPU Time of fewer DRTs (to target dimension 150) against size of datasets", "Dataset Size n", "CPU Execution Time (ns)", "img/drt_comparisons", X11 };

  /*
  //  L2, MaxDev and CPUTime for all datasets and all DRTs except APLA
  vector<vector<double>> dataset_samples;
  for (int di = 0; di<datasets.size(); di++) {
    dataset = parse_ucr_dataset(datasets[di], ucr_datasets_loc,  DatasetType::TRAIN_APPEND_TEST);
    if (dataset.size() < 10'000)
      continue;

    dataset.resize(10'000);
    for (int i=0; i<10; i++) {
      vector<double> dataset_sample(dataset.begin()+1000*i, dataset.begin()+1000*i+1000 );
      z_norm::z_normalise(dataset_sample);
      dataset_samples.push_back(dataset_sample);
    }

  }
  PlotDetails pd_l2 = { "Mean Euclidean Distance of DRTs on UCR datasets", "Target Dimension m", "L2"
    , "img/drt_comparisons/all_datasets/", PDF };
  PlotDetails pd_maxdev = { "Mean Maximum Deviation of DRTs on UCR datasets", "Target Dimension m", "Maximum Deviation"
    , "img/drt_comparisons/all_datasets/", PDF };

  PlotDetails pd_time = { "Mean CPU Time of DRTs on UCR datasets", "Target Dimension m", "CPU Execution Time (ns)"
    , "img/drt_comparisons/all_datasets/", PDF };
  //plot::plot_mean_bars_generated(dataset_samples, { 15, 30, 60, 90 }, l2_generators, pd_l2);
  //plot::plot_mean_bars_generated(dataset_samples, { 15, 30, 60, 90 }, maxdev_generators, pd_maxdev);
  plot::plot_mean_bars_generated(dataset_samples, { 15, 30, 60, 90 }, time_generators, pd_time);
  return 0;
  */


  /********************************** R Tree implementation **************************************/
  /*
  vector<unsigned int> divs = { 5, 13, 25, 109 };
  for (auto di : divs) {
    dataset = parse_ucr_dataset(datasets[di], ucr_datasets_loc,  DatasetType::TRAIN_APPEND_TEST);

    dataset.resize(20'000);
    z_norm::z_normalise(dataset);

    //z_norm::z_normalise(dataset);
    // choose to have 30 segments for subsequences of size 300
    const unsigned int seq_size = 600;
    const unsigned int NS1 = 15;
    const unsigned int NS2 = 30;
    const unsigned int NS3 = 60;

    vector<pla::APLA_DRT> apla_drts = { capla_eval::generate_mean_DRT_COMPR(5)
    , d_w_proj_apla_f_uncompr
    , rdp_f_uncompr
    , bottom_up_f_uncompr
    , sw_f_uncompr };
    for (int apla_i=0; apla_i < apla_drts.size(); apla_i++) {
      pla::APLA_DRT apla_drt = apla_drts[apla_i];
      double pp_ns1, pp_ns2, pp_ns3;

      auto retrieval_f = [](const unsigned int& i, const vector<double>& q) {
	return std::vector<std::array<const double*,2>>( {{ q.data()+i, q.data()+i+seq_size-1 }} );
      };
      // test on 100 trials
      unsigned int trial_incr = dataset.size() / 100;


      {

      RTree<apla_bounds::AplaMBR<NS1>, unsigned int> r_tree2(40,10 , apla_bounds::mbr_area<NS1> , apla_bounds::mbr_merge<NS1> , apla_bounds::dist_to_mbr_sqr<NS1>);

      //std::cout << "mbrs for apla NS1 calculating" << std::endl;
      auto vec2_of_mbrs = apla_bounds::vec_to_subseq_mbrs<NS1>(dataset,seq_size,apla_drt);

      //std::cout << "inserting" << std::endl;
      for (int i=0; i<vec2_of_mbrs.size(); i++) {
	r_tree2.insert( vec2_of_mbrs[i], i );
	if (i % (dataset.size()/10) == 0) {
	  //std::cout << (i / (dataset.size()/10)) << "0%..";
	}
      }
      //std::cout << std::endl;

      //std::cout << "searching" << std::endl;
      unsigned int trials = 0;
      for (int i=0;i<dataset.size()-seq_size; i+=trial_incr) {
	trials++;
	std::vector<double> query( dataset.begin()+i, dataset.begin()+seq_size+i);
	NormalFunctor noise(0,0.0,0.1);
	for (int i=0; i<seq_size; i++) { // perturb query a bit to make it new
	  query[i]+=noise();
	}
	pp_ns1 += r_tree2.pruning_power(query, retrieval_f, dataset);
      }
      //std::cout << "done NS1" << std::endl;
      pp_ns1/=trials;

      }

      {
      RTree<apla_bounds::AplaMBR<NS2>, unsigned int> r_tree2_NS2(40,10 , apla_bounds::mbr_area<NS2> , apla_bounds::mbr_merge<NS2> , apla_bounds::dist_to_mbr_sqr<NS2>);

      //std::cout << "mbrs for apla 11 calculating" << std::endl;
      auto vec2_of_mbrs_NS2 = apla_bounds::vec_to_subseq_mbrs<NS2>(dataset,seq_size,capla_eval::generate_mean_DRT_COMPR(5));

      //std::cout << "inserting" << std::endl;
      for (int i=0; i<vec2_of_mbrs_NS2.size(); i++) {
	r_tree2_NS2.insert( vec2_of_mbrs_NS2[i], i );
	if (i % (dataset.size()/10) == 0) {
	  //std::cout << (i / (dataset.size()/10)) << "0%..";
	}
      }
      //std::cout << std::endl;

      //std::cout << "searching" << std::endl;
      unsigned int trials = 0;
      for (int i=0;i<dataset.size()-seq_size; i+=trial_incr) {
	trials++;
	std::vector<double> query( dataset.begin()+i, dataset.begin()+seq_size+i);
	NormalFunctor noise(0,0.0,0.1);
	for (int i=0; i<seq_size; i++) { // perturb query a bit to make it new
	  query[i]+=noise();
	}

	pp_ns2 += r_tree2_NS2.pruning_power(query, retrieval_f, dataset);
      }
      //std::cout << "done NS2" << std::endl;
      pp_ns2/=trials;

      }

      {
      RTree<apla_bounds::AplaMBR<NS3>, unsigned int> r_tree2_NS3(40,10 , apla_bounds::mbr_area<NS3> , apla_bounds::mbr_merge<NS3> , apla_bounds::dist_to_mbr_sqr<NS3>);
      
      //std::cout << "mbrs for apla NS3 calculating" << std::endl;
      auto vec2_of_mbrs_NS3 = apla_bounds::vec_to_subseq_mbrs<NS3>(dataset,seq_size,capla_eval::generate_mean_DRT_COMPR(5));
      
      //std::cout << "inserting" << std::endl;
      for (int i=0; i<vec2_of_mbrs_NS3.size(); i++) {
	r_tree2_NS3.insert( vec2_of_mbrs_NS3[i], i );
	if (i % (dataset.size()/10) == 0) {
	  //std::cout << (i / (dataset.size()/10)) << "0%..";
	}
      }
      //std::cout << std::endl;

      //std::cout << "searching" << std::endl;
      unsigned int trials = 0;
      for (int i=0;i<dataset.size()-seq_size; i+=trial_incr) {
	trials++;
	std::vector<double> query( dataset.begin()+i, dataset.begin()+seq_size+i);
	NormalFunctor noise(0,0.0,0.1);
	for (int i=0; i<seq_size; i++) { // perturb query a bit to make it new
	  query[i]+=noise();
	}

	pp_ns3 += r_tree2_NS3.pruning_power(query, retrieval_f, dataset);
      }

      //std::cout << "done NS3" << std::endl;
      pp_ns3/=trials;
      
      //std::cout << "trials : "<< trials << std::endl;
      }

	std::cout << pp_ns1 << "," 
	<< pp_ns2 << "," 
	<< pp_ns3 << std::endl;
    }
  }
  */
  /********************************************************************************************/

  /************************* GEMINI vs SeqScan ************************************************/
  /*
  {
  RandomWalk walk( NormalFunctor(1) ); 
  walk.gen_steps(50'000);
  dataset = std::vector<double>(walk.get_walk().cbegin(), walk.get_walk().cend()); 

  std::cout << dataset.size() << std::endl;
  z_norm::z_normalise(dataset);


  const unsigned int seq_size = 600;
  const unsigned int NS1 = 30;
  RTree<apla_bounds::AplaMBR<NS1>, unsigned int> r_tree4(40,10 , apla_bounds::mbr_area<NS1> , apla_bounds::mbr_merge<NS1> , apla_bounds::dist_to_mbr_sqr<NS1>);

  auto retrieval_f = [](const unsigned int& i, const vector<double>& q) {
    return std::vector<std::array<const double*,2>>( {{ q.data()+i, q.data()+i+seq_size-1 }} );
  };

  std::cout << "calculating mbrs" << std::endl;
  auto vec4_of_mbrs = apla_bounds::vec_to_subseq_mbrs<NS1>(dataset,seq_size, rdp_f_uncompr);

  std::cout << "inserting" << std::endl;
  for (int i=0; i<vec4_of_mbrs.size(); i++) {
    r_tree4.insert( vec4_of_mbrs[i], i );
  }

  std::cout << "number of mbrs added " << vec4_of_mbrs.size() << std::endl;
  std::cout << "number of reachable leaves " << r_tree4.get_size_leaves() << std::endl;

  std::cout << "searching" << std::endl;
  vector<double> rtree_times(2);
  vector<double> seqscan_times(2);
  double trials = 0;
  for (int i=0;i<dataset.size()-seq_size; i+=10'000) {
    trials += 1;
    std::vector<double> query2( dataset.begin()+i, dataset.begin()+seq_size+i);
    std::vector<double> query( query2 );
    NormalFunctor noise(0,0.0,0.001);
    for (int i=0; i<seq_size; i++) { // perturb query a bit to make it new
      query[i]+=noise();
    }

    auto start_rtree = std::chrono::high_resolution_clock::now();
    r_tree4.sim_search(query, 0.1);
    auto end_rtree = std::chrono::high_resolution_clock::now();
    auto dur_rtree = std::chrono::duration_cast<std::chrono::microseconds>(end_rtree - start_rtree);
    rtree_times[1] += dur_rtree.count()/1000.0;

    start_rtree = std::chrono::high_resolution_clock::now();
    r_tree4.knn_search(query, 30, retrieval_f, dataset);
    end_rtree = std::chrono::high_resolution_clock::now();
    dur_rtree = std::chrono::duration_cast<std::chrono::microseconds>(end_rtree - start_rtree);
    rtree_times[0] += dur_rtree.count()/1000.0;

    auto start_seqscan = std::chrono::high_resolution_clock::now();
    seq_scan::find_similar_subseq_indexes(dataset, query, 0.1);
    auto end_seqscan = std::chrono::high_resolution_clock::now();
    auto dur_seqscan = std::chrono::duration_cast<std::chrono::microseconds>(end_seqscan - start_seqscan);
    seqscan_times[1] += dur_seqscan.count()/1000.0;

    start_seqscan = std::chrono::high_resolution_clock::now();
    seq_scan::find_k_closest_indexes(dataset, query, 30);
    end_seqscan = std::chrono::high_resolution_clock::now();
    dur_seqscan = std::chrono::duration_cast<std::chrono::microseconds>(end_seqscan - start_seqscan);
    seqscan_times[0] += dur_seqscan.count()/1000.0;

  }
  rtree_times[0] /= trials;
  rtree_times[1] /= trials;
  seqscan_times[0] /= trials;
  seqscan_times[1] /= trials;

  PlotDetails pd = { "Time taken by Sim Search and K-NN Methods", "", "Time (MS)", "img/sim_search/", PDF };
  Series rtree_series = {rtree_times, "RTree"};
  Series seqscan_series = {seqscan_times, "Sequential Scan"};
  vector<string> x_labels = { "K-NN", "Similarity Search" };
  vector<Series> vs = {rtree_series, seqscan_series};
  plot::barplot_many_series(vs, x_labels, pd);

  }
  */
  /*****************************************/

  /*
  di = 109; 
  dataset = parse_ucr_dataset(datasets[di], ucr_datasets_loc,  DatasetType::TRAIN_APPEND_TEST);

  dataset.resize(1000);
  z_norm::z_normalise(dataset);

  Seqd x;
  for (int i=0; i<dataset.size(); i++) {
    x.push_back(i);
  }

  Seqd x_subset;
  for (int i=150; i<231; i++)
    x_subset.push_back(i);
  Seqd data_subset(dataset.data() + 150, dataset.data() + 231);
  Line l1 = { x_subset, data_subset,  "original" };

  const unsigned int seq_size = 80;
  const unsigned int NS = 8;
  auto retrieval_f = [](const unsigned int& i, const vector<double>& q) {
    return std::vector<std::array<const double*,2>>( {{ q.data()+i, q.data()+i+seq_size-1 }} );
  };
  RTree<apla_bounds::AplaMBR<NS>, unsigned int> r_tree(40,10 , apla_bounds::mbr_area<NS> , apla_bounds::mbr_merge<NS> , apla_bounds::dist_to_mbr_sqr<NS>);
  auto vec_of_mbrs = apla_bounds::vec_to_subseq_mbrs<NS>(dataset,seq_size,bottom_up_f_uncompr);
  for (int i=0; i<vec_of_mbrs.size(); i++) {
    r_tree.insert(vec_of_mbrs[i], i);
  }
  Seqd query(dataset.begin() + 150, dataset.begin() + 231);
  Seqd dataset_x; for (int i=0; i<dataset.size(); i++) dataset_x.push_back(i);
  std::vector<Line> vl;
  Line dataset_line = { dataset_x, dataset, "ECG Data" };
  vl.push_back(dataset_line);
  auto sim_results = r_tree.knn_search(query, 15, retrieval_f, dataset);
  for (const auto& [l_ptr,r_ptr] : sim_results) {
    Seqd res(l_ptr, r_ptr+1);
    Seqd res_x; for (int i=0; i<res.size(); i++) res_x.push_back( l_ptr - dataset.data() + i );
    vl.push_back({res_x, res, "distance 2 away"});
  }
  vl.push_back(l1);
  PlotDetails p_sim = { "K-NN on Strawberry Dataset", "Time", "Value", "./img/sim_search/", PDF };
  plot::plot_lines(vl,p_sim);
  */


  /*
  const unsigned int seq_size = ;
  const unsigned int NS1 = 30;
  RTree<apla_bounds::AplaMBR<NS1>, unsigned int> r_tree4(40,10 , apla_bounds::mbr_area<NS1> , apla_bounds::mbr_merge<NS1> , apla_bounds::dist_to_mbr_sqr<NS1>);

  auto retrieval_f = [](const unsigned int& i, const vector<double>& q) {
    return std::vector<std::array<const double*,2>>( {{ q.data()+i, q.data()+i+seq_size-1 }} );
  };

  std::cout << "calculating mbrs" << std::endl;
  auto vec4_of_mbrs = apla_bounds::vec_to_subseq_mbrs<NS1>(dataset,seq_size, rdp_f_uncompr);

  std::cout << "inserting" << std::endl;
  for (int i=0; i<vec4_of_mbrs.size(); i++) {
    r_tree4.insert( vec4_of_mbrs[i], i );
  }
  */

  /********************/

  return 0;
}
