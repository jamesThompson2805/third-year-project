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
  unsigned int di = 13;
  vector<double> dataset = parse_ucr_dataset(datasets[di], ucr_datasets_loc,  DatasetType::TRAIN);
  std::cout << dataset.size() << std::endl;
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
  //     demo();
  
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
  vector<LineGenerator> l2_generators = {
	paa_gen_l2
	, pla_gen_l2
	, apca_gen_l2
	, d_w_apla_gen_l2
	, d_w_proj_apla_gen_l2
	, c_d_w_mean_s1_gen_l2
	, c_d_w_tri_gen_l2
	, apla_gen_l2
	//, rdp_gen_l2
	//, bot_gen_l2
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
  vector<LineGenerator> maxdev_generators = {
	paa_gen_maxdev
	, pla_gen_maxdev
	, apca_gen_maxdev
	, d_w_apla_gen_maxdev
	, d_w_proj_apla_gen_maxdev
	, c_d_w_mean_s1_gen_maxdev
	, c_d_w_tri_gen_maxdev
	, apla_gen_maxdev
	//, rdp_gen_maxdev
	//, bot_gen_maxdev
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
  vector<LineGenerator> time_generators = {
	paa_gen_time
	, apca_gen_time
	, pla_gen_time
	, d_w_apla_gen_time
	, d_w_proj_apla_gen_time
	, c_d_w_mean_s1_gen_time
	, c_d_w_tri_gen_time
	, apla_gen_time
	//, rdp_gen_time
	//, bot_gen_time
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
  LinePGGenerator sw_gen_time = { sliding_w_time, "Sliding Window" };
  LinePGGenerator swing_gen_time = { swing_time, "Swing" };

  vector<LinePGGenerator> compr_generators = { top_down_gen, bottom_up_gen, sw_gen, swing_gen };
  vector<LinePGGenerator> etime_generators = { top_down_gen_time, bottom_up_gen_time, sw_gen_time, swing_gen_time };

  vector<unsigned int> compr_ds = { 5, 13, 25, 109 };
  for (auto compr_i : compr_ds) {
    dataset = parse_ucr_dataset(datasets[compr_i], ucr_datasets_loc,  DatasetType::TRAIN_APPEND_TEST);
    std::cout << dataset.size() << std::endl;
    dataset.resize(20'000);
    z_norm::z_normalise(dataset);
    PlotDetails pd_compr = { "Compression Ratio of DRTs on " + datasets[compr_i] + " dataset", "Epsilon Value", "Compression Ratio", "img/drt_comparisons/", PDF };
    PlotDetails pd_etime = { "CPU Time of DRTs on " + datasets[compr_i] + " dataset", "Epsilon Value", "CPU Execution Time (ns)", "img/drt_comparisons/", PDF };
    plot::plot_barsPG_generated(dataset, { 0.05, 0.1, 0.15, 0.2, 0.25 }, compr_generators, pd_compr);
    plot::plot_barsPG_generated(dataset, { 0.05, 0.1, 0.15, 0.2, 0.25 }, etime_generators, pd_etime);
  }
  

  /**************/




  PlotDetails p_l2 = { "Euclidean Distance of DRTs on Chlorine dataset", "Target Dimension m", "L2", "img/drt_comparisons/", PDF };
  PlotDetails p_maxdev = { "Maximum Deviation of DRTs on Chlorine dataset", "Target Dimension m", "Maximum Deviation", "img/drt_comparisons/", PDF };
  PlotDetails p_time = { "CPU Time of DRTs on Synthetic dataset", "Target Dimension m", "CPU Execution Time (ns)", "img/drt_comparisons/", PDF };
  PlotDetails p_time_size = { "CPU Time of fewer DRTs (to target dimension 150) against size of datasets", "Dataset Size n", "CPU Execution Time (ns)", "img/drt_comparisons", X11 };

  /*
  vector<unsigned int> ds_indexes;
  for (int i=0; i<datasets.size(); i+=4) ds_indexes.push_back(i);
  vector<string> names; std::for_each(ds_indexes.begin(), ds_indexes.end(), [&](unsigned int i){ names.push_back(datasets[i]); });
  plot::plot_lines_generated_ucr_average(names, ucr_datasets_loc, 4096
      //, {300, 600, 900, 1200, 1500, 1800, 2100, 2400, 2700, 3000}
      , {30, 60, 90, 120, 150, 180, 210, 240, 270, 300}
      , l2_generators, p_l2);
  */

  //PlotDetails pd_l2 = { "Euclidean Distance of DRTs on Synthetic dataset", "Target Dimension m", "L2", "img/drt_comparisons/", PDF };
  //PlotDetails pd_maxdev = { "Maximum Deviation of DRTs on Synthetic dataset", "Target Dimension m", "Maximum Deviation", "img/drt_comparisons/", PDF };
  //PlotDetails pd_time = { "CPU Time of DRTs on Synthetic dataset", "Target Dimension m", "CPU Execution Time (ns)", "img/drt_comparisons/", PDF };
  //plot::plot_bars_generated(dataset2, { 60, 120, 180, 240 }, l2_generators, pd_l2);
  //plot::plot_bars_generated(dataset2, { 60, 120, 180, 240 }, maxdev_generators, pd_maxdev);
  //plot::plot_bars_generated(dataset2, { 60, 120, 180, 240 }, time_generators, pd_time);



  /********************************** R Tree implementation **************************************/
  // choose to have 30 segments for subsequences of size 300
  const unsigned int seq_size = 1024;
  const unsigned int NS = 21;
  const unsigned int NS2 = 11;
  const unsigned int NS3 = 5;

  double capla_pp_ns1, rdp_pp_ns1, b_u_pp_ns1;
  double capla_pp_ns2, rdp_pp_ns2, b_u_pp_ns2;
  double capla_pp_ns3, rdp_pp_ns3, b_u_pp_ns3;

  auto retrieval_f = [](const unsigned int& i, const vector<double>& q) {
    return std::vector<std::array<const double*,2>>( {{ q.data()+i, q.data()+i+seq_size-1 }} );
  };

  /*
  {

  RTree<apla_bounds::AplaMBR<NS>, unsigned int> r_tree(40,10 , apla_bounds::mbr_area<NS> , apla_bounds::mbr_merge<NS> , apla_bounds::dist_to_mbr_sqr<NS>);
  RTree<apla_bounds::AplaMBR<NS>, unsigned int> r_tree2(40,10 , apla_bounds::mbr_area<NS> , apla_bounds::mbr_merge<NS> , apla_bounds::dist_to_mbr_sqr<NS>);
  RTree<apla_bounds::AplaMBR<NS>, unsigned int> r_tree3(40,10 , apla_bounds::mbr_area<NS> , apla_bounds::mbr_merge<NS> , apla_bounds::dist_to_mbr_sqr<NS>);
  RTree<apla_bounds::AplaMBR<NS>, unsigned int> r_tree4(40,10 , apla_bounds::mbr_area<NS> , apla_bounds::mbr_merge<NS> , apla_bounds::dist_to_mbr_sqr<NS>);

  //auto vec_of_mbrs = apla_bounds::vec_to_subseq_mbrs<NS>(dataset,300,exact_dp::min_mse_pla);
  std::cout << "mbrs for capla 21 calculating" << std::endl;
  auto vec2_of_mbrs = apla_bounds::vec_to_subseq_mbrs<NS>(dataset,seq_size,capla_eval::generate_mean_DRT_COMPR(4));
  std::cout << "mbrs for rdp 21 calculating" << std::endl;
  auto vec3_of_mbrs = apla_bounds::vec_to_subseq_mbrs<NS>(dataset,seq_size,rdp_f_uncompr);
  std::cout << "mbrs for b-u 21 calculating" << std::endl;
  auto vec4_of_mbrs = apla_bounds::vec_to_subseq_mbrs<NS>(dataset,seq_size,bottom_up_f_uncompr);

  std::cout << "inserting" << std::endl;
  for (int i=0; i<vec2_of_mbrs.size(); i++) {
    //r_tree.insert( vec_of_mbrs[i], i );
    r_tree2.insert( vec2_of_mbrs[i], i );
    r_tree3.insert( vec3_of_mbrs[i], i );
    r_tree4.insert( vec4_of_mbrs[i], i );
  }

  std::cout << "searching" << std::endl;
  unsigned int trials = 0;
  for (int i=0;i<dataset.size()-seq_size; i+=741) {
    trials++;
    std::vector<double> query( dataset.begin()+i, dataset.begin()+seq_size+i);
    NormalFunctor noise(0,0.0,0.1);
    for (int i=0; i<seq_size; i++) { // perturb query a bit to make it new
      query[i]+=noise();
    }
    capla_pp_ns1 += r_tree2.pruning_power(query, retrieval_f, dataset);
    rdp_pp_ns1 += r_tree3.pruning_power(query, retrieval_f, dataset);
    b_u_pp_ns1 += r_tree4.pruning_power(query, retrieval_f, dataset);
  }
  std::cout << "done 21" << std::endl;
  capla_pp_ns1/=trials;
  rdp_pp_ns1/=trials;
  b_u_pp_ns1/=trials;

  }

  {
  RTree<apla_bounds::AplaMBR<NS2>, unsigned int> r_tree2_NS2(40,10 , apla_bounds::mbr_area<NS2> , apla_bounds::mbr_merge<NS2> , apla_bounds::dist_to_mbr_sqr<NS2>);
  RTree<apla_bounds::AplaMBR<NS2>, unsigned int> r_tree3_NS2(40,10 , apla_bounds::mbr_area<NS2> , apla_bounds::mbr_merge<NS2> , apla_bounds::dist_to_mbr_sqr<NS2>);
  RTree<apla_bounds::AplaMBR<NS2>, unsigned int> r_tree4_NS2(40,10 , apla_bounds::mbr_area<NS2> , apla_bounds::mbr_merge<NS2> , apla_bounds::dist_to_mbr_sqr<NS2>);

  std::cout << "mbrs for capla 11 calculating" << std::endl;
  auto vec2_of_mbrs_NS2 = apla_bounds::vec_to_subseq_mbrs<NS2>(dataset,seq_size,capla_eval::generate_mean_DRT_COMPR(4));
  std::cout << "mbrs for rdp 11 calculating" << std::endl;
  auto vec3_of_mbrs_NS2 = apla_bounds::vec_to_subseq_mbrs<NS2>(dataset,seq_size,rdp_f_uncompr);
  std::cout << "mbrs for b-u 11 calculating" << std::endl;
  auto vec4_of_mbrs_NS2 = apla_bounds::vec_to_subseq_mbrs<NS2>(dataset,seq_size,bottom_up_f_uncompr);

  std::cout << "inserting" << std::endl;
  for (int i=0; i<vec2_of_mbrs_NS2.size(); i++) {
    r_tree2_NS2.insert( vec2_of_mbrs_NS2[i], i );
    r_tree3_NS2.insert( vec3_of_mbrs_NS2[i], i );
    r_tree4_NS2.insert( vec4_of_mbrs_NS2[i], i );
  }

  std::cout << "searching" << std::endl;
  unsigned int trials = 0;
  for (int i=0;i<dataset.size()-seq_size; i+=741) {
    trials++;
    std::vector<double> query( dataset.begin()+i, dataset.begin()+seq_size+i);
    NormalFunctor noise(0,0.0,0.1);
    for (int i=0; i<seq_size; i++) { // perturb query a bit to make it new
      query[i]+=noise();
    }

    capla_pp_ns2 += r_tree2_NS2.pruning_power(query, retrieval_f, dataset);
    rdp_pp_ns2 += r_tree3_NS2.pruning_power(query, retrieval_f, dataset);
    b_u_pp_ns2 += r_tree4_NS2.pruning_power(query, retrieval_f, dataset);
  }
  std::cout << "done 11" << std::endl;
  capla_pp_ns2/=trials;
  rdp_pp_ns2/=trials;
  b_u_pp_ns2/=trials;

  }

  {
  RTree<apla_bounds::AplaMBR<NS3>, unsigned int> r_tree2_NS3(40,10 , apla_bounds::mbr_area<NS3> , apla_bounds::mbr_merge<NS3> , apla_bounds::dist_to_mbr_sqr<NS3>);
  RTree<apla_bounds::AplaMBR<NS3>, unsigned int> r_tree3_NS3(40,10 , apla_bounds::mbr_area<NS3> , apla_bounds::mbr_merge<NS3> , apla_bounds::dist_to_mbr_sqr<NS3>);
  RTree<apla_bounds::AplaMBR<NS3>, unsigned int> r_tree4_NS3(40,10 , apla_bounds::mbr_area<NS3> , apla_bounds::mbr_merge<NS3> , apla_bounds::dist_to_mbr_sqr<NS3>);
  
  std::cout << "mbrs for capla 5 calculating" << std::endl;
  auto vec2_of_mbrs_NS3 = apla_bounds::vec_to_subseq_mbrs<NS3>(dataset,seq_size,capla_eval::generate_mean_DRT_COMPR(4));
  std::cout << "mbrs for rdp 5 calculating" << std::endl;
  auto vec3_of_mbrs_NS3 = apla_bounds::vec_to_subseq_mbrs<NS3>(dataset,seq_size,rdp_f_uncompr);
  std::cout << "mbrs for b-u 5 calculating" << std::endl;
  auto vec4_of_mbrs_NS3 = apla_bounds::vec_to_subseq_mbrs<NS3>(dataset,seq_size,bottom_up_f_uncompr);
  
  std::cout << "inserting" << std::endl;
  for (int i=0; i<vec2_of_mbrs_NS3.size(); i++) {
    r_tree2_NS3.insert( vec2_of_mbrs_NS3[i], i );
    r_tree3_NS3.insert( vec3_of_mbrs_NS3[i], i );
    r_tree4_NS3.insert( vec4_of_mbrs_NS3[i], i );
  }

  std::cout << "searching" << std::endl;
  unsigned int trials = 0;
  for (int i=0;i<dataset.size()-seq_size; i+=741) {
    trials++;
    std::vector<double> query( dataset.begin()+i, dataset.begin()+seq_size+i);
    NormalFunctor noise(0,0.0,0.1);
    for (int i=0; i<seq_size; i++) { // perturb query a bit to make it new
      query[i]+=noise();
    }

    capla_pp_ns3 += r_tree2_NS3.pruning_power(query, retrieval_f, dataset);
    rdp_pp_ns3 += r_tree3_NS3.pruning_power(query, retrieval_f, dataset);
    b_u_pp_ns3 += r_tree4_NS3.pruning_power(query, retrieval_f, dataset);
  }

  std::cout << "done 5" << std::endl;
  capla_pp_ns3/=trials;
  rdp_pp_ns3/=trials;
  b_u_pp_ns3/=trials;

  }

  std::cout << "capla pruning powers " << NS3 << " : " << capla_pp_ns3 << ", " << NS2 << " : " << capla_pp_ns2 << ", " << NS << " : " << capla_pp_ns1 << std::endl;
  std::cout << "rdp pruning powers " << NS3 << " : " << rdp_pp_ns3 << ", " << NS2 << " : " << rdp_pp_ns2 << ", " << NS << " : " << rdp_pp_ns1 << std::endl;
  std::cout << "b-u pruning powers " << NS3 << " : " << b_u_pp_ns3 << ", " << NS2 << " : " << b_u_pp_ns2 << ", " << NS << " : " << b_u_pp_ns1 << std::endl;
  */

  /*
  Series capla_pruning = { { 0.0717536, 0.000605522, 2.71927e-05 }, "Sliding Window" };
  Series rdp_pruning = { { 0.00894341, 5.63206e-05, 2.86695e-05 }, "RDP" };
  Series b_u_pruning = { { 0.00896001, 5.49966e-05, 2.87714e-05 }, "B-U" };
  vector<Series> vs = { capla_pruning, rdp_pruning, b_u_pruning };
  PlotDetails pd = { "Pruning Power of PLA Methods", "", "Pruning Power", "img/", PDF };
  plot::barplot_many_series(vs, pd);
  */


  {/*
  RTree<apla_bounds::AplaMBR<NS>, unsigned int> r_tree4(40,10 , apla_bounds::mbr_area<NS> , apla_bounds::mbr_merge<NS> , apla_bounds::dist_to_mbr_sqr<NS>);

  std::cout << "calculating mbrs" << std::endl;
  auto vec4_of_mbrs = apla_bounds::vec_to_subseq_mbrs<NS>(dataset,seq_size,bottom_up_f_uncompr);

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

  PlotDetails pd = { "Time taken by Sim Search and K-NN Methods", "", "Time (MS)", "img/", PDF };
  Series rtree_series = {rtree_times, "RTree"};
  Series seqscan_series = {seqscan_times, "Sequential Scan"};
  vector<Series> vs = {rtree_series, seqscan_series};
  plot::barplot_many_series(vs, pd);

  */
  }

  /********************/

  dataset2.resize(120);
  z_norm::z_normalise(dataset2);
  CauchyFunctor noise(7, 0.0, 1.0);
  vector<double> dataset3 = {0.0};
  for (int i=1; i<dataset2.size(); i++) {
    dataset3.push_back( dataset3[i-1] + noise());
  }
  z_norm::z_normalise(dataset3);
  Series s2 = { dataset3, "Standard Cauchy Walk" };
  //

  return 0;
}
