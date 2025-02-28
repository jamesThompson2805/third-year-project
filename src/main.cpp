#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <tuple>
#include <algorithm>

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

#include "conv_double_window.h"

#include "plotting/series_plotting.h"
#include "plotting/plot_dimreduct_paa.cpp"
#include "plotting/plot_dimreduct_pla.cpp"

#include "evaluations/general.h"
#include "evaluations/capla.h"

#include "random_walk.h"

#include "r_tree.h"
#include "lower_bounds_apla.h"


using std::vector, std::string, std::tuple;

int main()
{
  using namespace ucr_parsing;

  string ucr_datasets_loc = "external/data/UCRArchive_2018/";
  vector<string> datasets = parse_folder_names(ucr_datasets_loc);

  /**
  for (int i=0; i<datasets.size(); i++) {
    vector<double> dataset = parse_ucr_dataset(datasets[i], ucr_datasets_loc,  DatasetType::TRAIN);
    std::cout << i << " : " << datasets[i] << " : " << dataset.size() << std::endl;
  }
  **/

  // Favourites : 5 is Arrowhead, 13 is Chlorine, 28 is ECG200, 39 is Fifty Words, 128 is yoga, 46 is GuestureMidAirD1
  unsigned int di = 39;
  vector<double> dataset = parse_ucr_dataset(datasets[di], ucr_datasets_loc,  DatasetType::TRAIN);
  dataset.resize(400);
  z_norm::z_normalise(dataset);

  auto paa_f = [](const vector<double>& s, unsigned int num_params){ return paa::paa_to_seq( paa::paa(s, num_params), (s.size() / num_params) + (s.size() % num_params != 0)); };
  auto pla_f = [](const vector<double>& s, unsigned int num_params){ return pla::pla_to_seq(pla::pla(s, num_params), (2*s.size() / num_params) + (2*s.size() % num_params != 0)); };

  auto d_w_apca_f = [](const vector<double>& s, unsigned int num_params){ return d_w::simple_paa(s, num_params, 5, 5); };
  auto d_w_proj_apca_f = [](const vector<double>& s, unsigned int num_params){ return d_w::y_proj_paa(s, num_params, 5, 5); };

  auto d_w_apla_f = [](const vector<double>& s, unsigned int num_params){ return pla::apla_to_seq(d_w::simple_pla(s, num_params, 3, 3)); };
  auto d_w_proj_apla_f = [](const vector<double>& s, unsigned int num_params){ return pla::apla_to_seq(d_w::y_proj_pla(s, num_params, 5, 5)); };

  auto exact_apaa_f = [](const vector<double>& s, unsigned int num_params){ return paa::apca_to_seq(exact_dp::min_mse_paa(s, num_params)); }; 
  auto exact_apla_f = [](const vector<double>& s, unsigned int num_params){ return pla::apla_to_seq(exact_dp::min_mse_pla(s, num_params)); }; 

  auto apca_f = [](const vector<double>& s, unsigned int num_params){ return paa::apca_to_seq(apca::apca(s, num_params)); };

  auto rdp_f = [&](const Seqd& s, unsigned int parameter){ 
    auto apla = dac_curve_fitting::dac_linear(s, 0.1);
    if (apla.size() < parameter/3) segmerge::segment_to_dim(s,apla,parameter);
    if (apla.size() > parameter/3) segmerge::merge_to_dim(s,apla,parameter);
    return pla::apla_to_seq(apla);
  };
  auto bottom_up_f = [&](const Seqd& s, unsigned int parameter){ 
    auto apla = bottom_up::bottom_up(s, 0.1, bottom_up::se);
    if (apla.size() < parameter/3) segmerge::segment_to_dim(s,apla,parameter);
    if (apla.size() > parameter/3) segmerge::merge_to_dim(s,apla,parameter);
    return pla::apla_to_seq(apla);
  };

  vector<double> ldist = { 1.0/3.0, 1.0/3.0, 1.0/3.0};
  vector<double> rdist = { 0.0, 1.0/2.0, 1.0/2.0};
  auto conv_apla_f = [&ldist, &rdist](const vector<double>& s, unsigned int num_params){ return c_d_w::conv_pla(s, num_params, ldist, rdist); }; 

  RandomWalk walk( NormalFunctor(1) ); 
  walk.gen_steps(120);
  walk.save_walk("./tsv/testwalk1.tsv");
  vector<double> dataset2 = parse_tsv("./tsv/testwalk1.tsv",-1);
  z_norm::z_normalise(dataset2);
  Series s = { dataset2, "Normal Walk" };
  // plot_series(s, "./img/");


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

  /************************** MSE Evaluation against parameters *****************************/
  // PAA
  auto paa_gen_mse_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s, parameter, paa_f); };
  LineGenerator paa_gen_mse = { paa_gen_mse_f, "PAA" };
  // APCA
  auto apca_gen_mse_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s, parameter, apca_f); };
  LineGenerator apca_gen_mse = { apca_gen_mse_f, "APCA" };
  // PLA
  auto pla_gen_mse_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s, parameter, pla_f); };
  LineGenerator pla_gen_mse = { pla_gen_mse_f, "PLA" };
  // Double Window PLA and Proj
  auto d_w_apla_gen_mse_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s, parameter, d_w_apla_f); };
  LineGenerator d_w_apla_gen_mse = { d_w_apla_gen_mse_f, "Double Window PLA" };
  auto d_w_proj_apla_gen_mse_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s, parameter, d_w_proj_apla_f); };
  LineGenerator d_w_proj_apla_gen_mse = { d_w_proj_apla_gen_mse_f, "Projection PLA" };
  // CAPLA hypotheses
  auto c_d_w_mean_gen_mse_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s, parameter, capla_mean_f); };
  LineGenerator c_d_w_mean_gen_mse = { c_d_w_mean_gen_mse_f, "Faster Double Window PLA" };
  auto c_d_w_tri_gen_mse_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s, parameter, capla_tri_f); };
  LineGenerator c_d_w_tri_gen_mse = { c_d_w_tri_gen_mse_f, "Triangular Double Window PLA" };
  auto c_d_w_mean_s1_gen_mse_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s, parameter, capla_mean_s1_f); };
  LineGenerator c_d_w_mean_s1_gen_mse = { c_d_w_mean_s1_gen_mse_f, "mean s1 double window pla" };
  auto c_d_w_tri_s1_gen_mse_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s, parameter, capla_tri_s1_f); };
  LineGenerator c_d_w_tri_s1_gen_mse = { c_d_w_tri_s1_gen_mse_f, "Tri s1 double window pla" };
  auto apla_gen_mse_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s, parameter, exact_apla_f); };
  LineGenerator apla_gen_mse = { apla_gen_mse_f, "DP Exact PLA" };
  auto rdp_gen_mse_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s,parameter,rdp_f); };
  LineGenerator rdp_gen_mse = { rdp_gen_mse_f, "RDP" };
  auto bot_gen_mse_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s,parameter,bottom_up_f); };
  LineGenerator bot_gen_mse = { bot_gen_mse_f, "B-U" };
  vector<LineGenerator> mse_generators = {
	paa_gen_mse
	, apca_gen_mse
	, pla_gen_mse
	, d_w_apla_gen_mse
	, c_d_w_mean_gen_mse
	, c_d_w_tri_gen_mse
	, apla_gen_mse
	, rdp_gen_mse
	, bot_gen_mse
  };
  /**************/

  /************************** Euclidean Evaluation against parameters *****************************/
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
  LineGenerator d_w_apla_gen_l2 = { d_w_apla_gen_l2_f, "Double Window PLA" };
  auto d_w_proj_apla_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s, parameter, d_w_proj_apla_f); };
  LineGenerator d_w_proj_apla_gen_l2 = { d_w_proj_apla_gen_l2_f, "Projection PLA" };
  // CAPLA hypotheses
  auto c_d_w_mean_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s, parameter, capla_mean_f); };
  LineGenerator c_d_w_mean_gen_l2 = { c_d_w_mean_gen_l2_f, "Faster Double Window PLA" };
  auto c_d_w_tri_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s, parameter, capla_tri_f); };
  LineGenerator c_d_w_tri_gen_l2 = { c_d_w_tri_gen_l2_f, "Triangular Double Window PLA" };
  auto c_d_w_mean_s1_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s, parameter, capla_mean_s1_f); };
  LineGenerator c_d_w_mean_s1_gen_l2 = { c_d_w_mean_s1_gen_l2_f, "mean s1 double window pla" };
  auto c_d_w_tri_s1_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s, parameter, capla_tri_s1_f); };
  LineGenerator c_d_w_tri_s1_gen_l2 = { c_d_w_tri_s1_gen_l2_f, "Tri s1 double window pla" };
  auto apla_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s, parameter, exact_apla_f); };
  LineGenerator apla_gen_l2 = { apla_gen_l2_f, "DP Exact PLA" };
  auto rdp_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s,parameter,rdp_f); };
  LineGenerator rdp_gen_l2 = { rdp_gen_l2_f, "RDP" };
  auto bot_gen_l2_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::l2_of_method(s,parameter,bottom_up_f); };
  LineGenerator bot_gen_l2 = { bot_gen_l2_f, "B-U" };
  vector<LineGenerator> l2_generators = {
	paa_gen_l2
	, apca_gen_l2
	, pla_gen_l2
	, d_w_apla_gen_l2
	, c_d_w_mean_gen_l2
	, c_d_w_tri_gen_l2
	, apla_gen_l2
	, rdp_gen_l2
	, bot_gen_l2
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
  LineGenerator d_w_apla_gen_maxdev = { d_w_apla_gen_maxdev_f, "Double Window PLA" };
  auto d_w_proj_apla_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, d_w_proj_apla_f); };
  LineGenerator d_w_proj_apla_gen_maxdev = { d_w_proj_apla_gen_maxdev_f, "Projection PLA" };
  // CAPLA hypotheses
  auto c_d_w_mean_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, capla_mean_f); };
  LineGenerator c_d_w_mean_gen_maxdev = { c_d_w_mean_gen_maxdev_f, "Faster Double Window PLA" };
  auto c_d_w_tri_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, capla_tri_f); };
  LineGenerator c_d_w_tri_gen_maxdev = { c_d_w_tri_gen_maxdev_f, "Triangular Double Window PLA" };
  auto c_d_w_mean_s1_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, capla_mean_s1_f); };
  LineGenerator c_d_w_mean_s1_gen_maxdev = { c_d_w_mean_s1_gen_maxdev_f, "mean s1 double window pla" };
  auto c_d_w_tri_s1_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, capla_tri_s1_f); };
  LineGenerator c_d_w_tri_s1_gen_maxdev = { c_d_w_tri_s1_gen_maxdev_f, "Tri s1 double window pla" };
  auto apla_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, exact_apla_f); };
  LineGenerator apla_gen_maxdev = { apla_gen_maxdev_f, "DP Exact PLA" };
  auto rdp_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, rdp_f); };
  LineGenerator rdp_gen_maxdev = { rdp_gen_maxdev_f, "RDP" };
  auto bot_gen_maxdev_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s,parameter,bottom_up_f); };
  LineGenerator bot_gen_maxdev = { bot_gen_maxdev_f, "B-U" };
  vector<LineGenerator> maxdev_generators = {
	paa_gen_maxdev
	, apca_gen_maxdev
	, pla_gen_maxdev
	, d_w_apla_gen_maxdev
	, c_d_w_mean_gen_maxdev
	, c_d_w_tri_gen_maxdev
	//, apla_gen_maxdev
	, rdp_gen_maxdev
	, bot_gen_maxdev
  };
  /**************/

  /************************** Time Evaluation against size *****************************/
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
  LineGenerator d_w_apla_gen_time = { d_w_apla_gen_time_f, "Double Window PLA" };
  auto d_w_proj_apla_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, d_w_proj_apla_f); };
  LineGenerator d_w_proj_apla_gen_time = { d_w_proj_apla_gen_time_f, "Projection PLA" };
  // CAPLA hypotheses
  auto c_d_w_mean_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, capla_mean_f); };
  LineGenerator c_d_w_mean_gen_time = { c_d_w_mean_gen_time_f, "Faster Double Window PLA" };
  auto c_d_w_tri_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, capla_tri_f); };
  LineGenerator c_d_w_tri_gen_time = { c_d_w_tri_gen_time_f, "Triangular Double Window PLA" };
  auto c_d_w_mean_s1_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, capla_mean_s1_f); };
  LineGenerator c_d_w_mean_s1_gen_time = { c_d_w_mean_s1_gen_time_f, "mean s1 double window pla" };
  auto c_d_w_tri_s1_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, capla_tri_s1_f); };
  LineGenerator c_d_w_tri_s1_gen_time = { c_d_w_tri_s1_gen_time_f, "Tri s1 double window pla" };
  auto apla_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, exact_apla_f); };
  LineGenerator apla_gen_time = { apla_gen_time_f, "DP Exact PLA" };
  auto rdp_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, rdp_f); };
  LineGenerator rdp_gen_time = { rdp_gen_time_f, "RDP" };
  auto bot_gen_time_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s,parameter,bottom_up_f); };
  LineGenerator bot_gen_time = { bot_gen_time_f, "B-U" };
  vector<LineGenerator> time_generators = {
	paa_gen_time
	, apca_gen_time
	, pla_gen_time
	, d_w_apla_gen_time
	, c_d_w_mean_gen_time
	, c_d_w_tri_gen_time
	, apla_gen_time
	, rdp_gen_time
	, bot_gen_time
  };
  /**************/

  /************************** time Evaluation against parameters *****************************/
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
    return general_eval::cputime_ms_of_method(s1, tdim,bottom_up_f);
  };
  LineGenerator bot_gen_size_time = { bot_gen_size_time_f, "B-U" };
  vector<LineGenerator> time_size_generators = {
	paa_gen_size_time
	, apca_gen_size_time
	, pla_gen_size_time
	, d_w_apla_gen_size_time
	, c_d_w_mean_gen_size_time
	, c_d_w_tri_gen_size_time
	//, apla_gen_size_time
	//, rdp_gen_size_time
	//, bot_gen_size_time
  };
  /**************/

  PlotDetails p_mse = { "MSE of DRTs against compressed dimension", "Target Dimension m", "MSE", "img/", X11 };
  PlotDetails p_l2 = { "Euclidean Distance of DRTs against compressed dimension", "Target Dimension m", "L2", "img/", X11 };
  PlotDetails p_maxdev = { "Maximum Deviation of DRTs against compressed dimension", "Target Dimension m", "Maximum Deviation", "img/", X11 };
  PlotDetails p_time = { "CPU Time of DRTs against compressed dimension", "Target Dimension m", "CPU Execution Time (ns)", "img/", X11 };
  PlotDetails p_time_size = { "CPU Time of DRTs against size of dataset", "Dataset Size n", "CPU Execution Time (ns)", "img/", X11 };
  vector<unsigned int> ds_indexes = { 5, 13, 28, 39, 46, 128 };
  vector<string> names; std::for_each(ds_indexes.begin(), ds_indexes.end(), [&](unsigned int i){ names.push_back(datasets[i]); });
  plot::plot_lines_generated_ucr_average(names, ucr_datasets_loc, 300
      , {300, 600, 900, 1200, 1500, 1800, 2100, 2400, 2700, 3000}
      , time_size_generators, p_time_size);

  /********************************** R Tree implementation **************************************
  // choose to have 3 segments for subsequences of size 30
  RTree<apla_bounds::AplaMBR<3>, unsigned int> r_tree(40,10
      , apla_bounds::mbr_area<3>
      , apla_bounds::mbr_merge<3>
      , apla_bounds::dist_to_mbr_sqr<3>);
  
  auto vec_of_mbrs = apla_bounds::vec_to_subseq_mbrs<3>(dataset,30,exact_dp::min_mse_pla);
  
  for (int i=0; i<vec_of_mbrs.size(); i++) {
    std::vector<double> query( dataset.begin()+i, dataset.begin()+30+i);
    if ( apla_bounds::dist_to_mbr_sqr<3>(query, vec_of_mbrs[i]) >= 1e-30) {
      std::cout << " i " << i << " wrong dist " << apla_bounds::dist_to_mbr_sqr<3>(query, vec_of_mbrs[i]) << std::endl;
    }

    r_tree.insert( vec_of_mbrs[i], i );
  }

  RTree<apla_bounds::AplaMBR<3>, unsigned int> r_tree2(40,10
      , apla_bounds::mbr_area<3>
      , apla_bounds::mbr_merge<3>
      , apla_bounds::dist_to_mbr_sqr<3>);
  
  auto vec2_of_mbrs = apla_bounds::vec_to_subseq_mbrs<3>(dataset,30,capla_eval::generate_mean_DRT_COMPR(4));
  
  for (int i=0; i<vec2_of_mbrs.size(); i++) {
    std::vector<double> query( dataset.begin()+i, dataset.begin()+30+i);
    if ( apla_bounds::dist_to_mbr_sqr<3>(query, vec2_of_mbrs[i]) >= 1e-30) {
      std::cout << " i " << i << " wrong dist " << apla_bounds::dist_to_mbr_sqr<3>(query, vec2_of_mbrs[i]) << std::endl;
    }

    r_tree2.insert( vec2_of_mbrs[i], i );
  }
  

  unsigned int t_ind = 0;
  std::vector<double> query0( dataset.begin()+t_ind, dataset.begin()+30+t_ind);
  for (int i=0; i<query0.size(); i++) { 
    std::cout << i << ": " << query0[i] << "	";
    if (i % 5 == 0)
      std::cout << "\n";
  }
  std::cout << "\n" <<  std::endl;
  std::cout << apla_bounds::dist_to_mbr_sqr<3>(query0, apla_bounds::vec_to_mbr<3>(query0, capla_eval::generate_mean_DRT_COMPR(4)) ) << std::endl;

  auto retrieval_f = [](const unsigned int& i, const vector<double>& q) {
    return std::vector<std::array<const double*,2>>( {{ q.data()+i, q.data()+i+29 }} );
  };
  std::vector<std::array<const double*,2>> k_nn = r_tree.knn_search(query0, 10, retrieval_f, dataset);
  std::cout << "best pruning power " << r_tree.pruning_power(query0, retrieval_f, dataset) << std::endl;
  for (auto& [ptr1,ptr2] : k_nn) {
    int diff = ptr1 - dataset.data();
    std::cout << "indexes : " << diff << std::endl;
  }
  std::cout << std::endl;
  std::cout << "capla pruning power " << r_tree2.pruning_power(query0, retrieval_f, dataset) << std::endl;
  std::vector<std::array<const double*,2>> k_n2 = r_tree.knn_search(query0, 10, retrieval_f, dataset);
  for (auto& [ptr1,ptr2] : k_n2) {
    int diff = ptr1 - dataset.data();
    std::cout << "indexes : " << diff << std::endl;
  }
  ********************/


  

  /*
  vector<double> mse_s;
  for (int j=4; j<=32; j*=2) {
    vector<double> ldist_test;
    vector<double> rdist_test;
    for (int i=0; i<j; i++) {
      ldist_test.push_back( 1.0 / (double) j );
      rdist_test.push_back( 1.0 / (double) j );
    }
    int datasets_assessed = 0;
    double avg_mse = 0.0;
    for (int i=0; i<datasets.size(); i++) {
      vector<double> dataset = parse_ucr_dataset(datasets[i], ucr_datasets_loc,  DatasetType::TEST_APPEND_TRAIN);
      z_norm::z_normalise(dataset);
      if (dataset.size() < 10'000) continue;
      vector<double> capla = pla::apla_to_seq( c_d_w::conv_pla(dataset, 120, ldist_test, rdist_test) );
      avg_mse += mse::se_between_seq(dataset, capla);
      datasets_assessed++;
    }
    avg_mse /= datasets_assessed;
    mse_s.push_back(avg_mse);
  }
  Series sa = { mse_s, "MSE for DW APLA" };
  plot_series(sa, "img/");
  */

  /*
  NormalFunctor noise(5, 0.0, 1.0);
  vector<double> dataset3;
  for (double d : dataset2) {
    dataset3.push_back( d + 2.0 + noise());
  }
  Series s2 = { dataset3, "Altered Walk" };
  //plot_series_diff( s, s2, "./img/");
  //
  */

  return 0;
}
