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
#include "exact_dp.h"

#include "conv_double_window.h"

#include "plotting/series_plotting.h"

#include "evaluations/general.h"
#include "evaluations/capla.h"

#include "random_walk.h"

using std::vector, std::string, std::tuple;

int main()
{
  using namespace ucr_parsing;

  string ucr_datasets_loc = "external/data/UCRArchive_2018/";
  vector<string> datasets = parse_folder_names(ucr_datasets_loc);

  /*
  for (int i=0; i<datasets.size(); i++) {
    vector<double> dataset = parse_ucr_dataset(datasets[i], ucr_datasets_loc,  DatasetType::TRAIN);
    z_norm::z_normalise(dataset);
    std::cout << datasets[i] << " : " << dataset.size() << std::endl;
  }
  */

  unsigned int di = 26;
  vector<double> dataset = parse_ucr_dataset(datasets[di], ucr_datasets_loc,  DatasetType::TRAIN);
  dataset.resize(40000);
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

  vector<double> ldist = { 1.0/3.0, 1.0/3.0, 1.0/3.0};
  vector<double> rdist = { 0.0, 1.0/2.0, 1.0/2.0};
  auto conv_apla_f = [&ldist, &rdist](const vector<double>& s, unsigned int num_params){ return c_d_w::conv_pla(s, num_params, ldist, rdist); }; 

  /************************** MSE Evaluation against parameters *****************************/
  { 
  PlotDetails p = { "title", "number of parameters", "mse", "no file path", X11 };
  // PAA
  auto paa_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s, parameter, paa_f); };
  LineGenerator paa_gen = { paa_gen_f, "paa" };
  // APCA
  auto apca_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s, parameter, apca_f); };
  LineGenerator apca_gen = { apca_gen_f, "apca" };
  // PLA
  auto pla_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s, parameter, pla_f); };
  LineGenerator pla_gen = { pla_gen_f, "pla" };
  // Double Window PLA and Proj
  auto d_w_apla_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s, parameter, d_w_apla_f); };
  LineGenerator d_w_apla_gen = { d_w_apla_gen_f, "double window pla" };
  auto d_w_proj_apla_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s, parameter, d_w_proj_apla_f); };
  LineGenerator d_w_proj_apla_gen = { d_w_proj_apla_gen_f, "projection pla" };
  // CAPLA hypotheses
  auto capla_mean_f =  capla_eval::generate_mean_DRT(5);
  auto c_d_w_mean_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s, parameter, capla_mean_f); };
  LineGenerator c_d_w_mean_gen = { c_d_w_mean_gen_f, "mean double window pla" };
  auto capla_tri_f =  capla_eval::generate_tri_DRT(5);
  auto c_d_w_tri_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s, parameter, capla_tri_f); };
  LineGenerator c_d_w_tri_gen = { c_d_w_tri_gen_f, "Tri double window pla" };
  auto capla_mean_s1_f =  capla_eval::generate_mean_skip_one_DRT(5);
  auto c_d_w_mean_s1_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s, parameter, capla_mean_s1_f); };
  LineGenerator c_d_w_mean_s1_gen = { c_d_w_mean_s1_gen_f, "mean s1 double window pla" };
  auto capla_tri_s1_f =  capla_eval::generate_tri_skip_one_DRT(5);
  auto c_d_w_tri_s1_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::mse_of_method(s, parameter, capla_tri_s1_f); };
  LineGenerator c_d_w_tri_s1_gen = { c_d_w_tri_s1_gen_f, "Tri s1 double window pla" };
  plot::plot_lines_generated(dataset
      , {15, 30, 45, 60, 75, 90}
      , {
	paa_gen
	, apca_gen
	, pla_gen
	, d_w_apla_gen
	, d_w_proj_apla_gen
	, c_d_w_mean_gen
	, c_d_w_mean_s1_gen
	, c_d_w_tri_gen
	, c_d_w_tri_s1_gen
      }, p);
  }

  /************************** maxdev Evaluation against parameters *****************************/
  {
  PlotDetails p = { "title", "number of parameters", "maxdev", "no file path", X11 };
  // PAA
  auto paa_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, paa_f); };
  LineGenerator paa_gen = { paa_gen_f, "paa" };
  // APCA
  auto apca_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, apca_f); };
  LineGenerator apca_gen = { apca_gen_f, "apca" };
  // PLA
  auto pla_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, pla_f); };
  LineGenerator pla_gen = { pla_gen_f, "pla" };
  // Double Window PLA and Proj
  auto d_w_apla_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, d_w_apla_f); };
  LineGenerator d_w_apla_gen = { d_w_apla_gen_f, "double window pla" };
  auto d_w_proj_apla_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, d_w_proj_apla_f); };
  LineGenerator d_w_proj_apla_gen = { d_w_proj_apla_gen_f, "projection pla" };
  // CAPLA hypotheses
  auto capla_mean_f =  capla_eval::generate_mean_DRT(5);
  auto c_d_w_mean_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, capla_mean_f); };
  LineGenerator c_d_w_mean_gen = { c_d_w_mean_gen_f, "mean double window pla" };
  auto capla_tri_f =  capla_eval::generate_tri_DRT(5);
  auto c_d_w_tri_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, capla_tri_f); };
  LineGenerator c_d_w_tri_gen = { c_d_w_tri_gen_f, "Tri double window pla" };
  auto capla_mean_s1_f =  capla_eval::generate_mean_skip_one_DRT(5);
  auto c_d_w_mean_s1_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, capla_mean_s1_f); };
  LineGenerator c_d_w_mean_s1_gen = { c_d_w_mean_s1_gen_f, "mean s1 double window pla" };
  auto capla_tri_s1_f =  capla_eval::generate_tri_skip_one_DRT(5);
  auto c_d_w_tri_s1_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::maxdev_of_method(s, parameter, capla_tri_s1_f); };
  LineGenerator c_d_w_tri_s1_gen = { c_d_w_tri_s1_gen_f, "Tri s1 double window pla" };
  plot::plot_lines_generated(dataset
      , {15, 30, 45, 60, 75, 90}
      , {
	paa_gen
	, apca_gen
	, pla_gen
	, d_w_apla_gen
	, d_w_proj_apla_gen
	, c_d_w_mean_gen
	, c_d_w_mean_s1_gen
	, c_d_w_tri_gen
	, c_d_w_tri_s1_gen
      }, p);
  }
  
  /************************** time Evaluation against parameters *****************************/
  {
  PlotDetails p = { "title", "number of parameters", "time", "no file path", X11 };
  // PAA
  auto paa_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, paa_f); };
  LineGenerator paa_gen = { paa_gen_f, "paa" };
  // APCA
  auto apca_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, apca_f); };
  LineGenerator apca_gen = { apca_gen_f, "apca" };
  // PLA
  auto pla_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, pla_f); };
  LineGenerator pla_gen = { pla_gen_f, "pla" };
  // Double Window PLA and Proj
  auto d_w_apla_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, d_w_apla_f); };
  LineGenerator d_w_apla_gen = { d_w_apla_gen_f, "double window pla" };
  auto d_w_proj_apla_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, d_w_proj_apla_f); };
  LineGenerator d_w_proj_apla_gen = { d_w_proj_apla_gen_f, "projection pla" };
  // CAPLA hypotheses
  auto capla_mean_f =  capla_eval::generate_mean_DRT(5);
  auto c_d_w_mean_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, capla_mean_f); };
  LineGenerator c_d_w_mean_gen = { c_d_w_mean_gen_f, "mean double window pla" };
  auto capla_tri_f =  capla_eval::generate_tri_DRT(5);
  auto c_d_w_tri_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, capla_tri_f); };
  LineGenerator c_d_w_tri_gen = { c_d_w_tri_gen_f, "Tri double window pla" };
  auto capla_mean_s1_f =  capla_eval::generate_mean_skip_one_DRT(5);
  auto c_d_w_mean_s1_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, capla_mean_s1_f); };
  LineGenerator c_d_w_mean_s1_gen = { c_d_w_mean_s1_gen_f, "mean s1 double window pla" };
  auto capla_tri_s1_f =  capla_eval::generate_tri_skip_one_DRT(5);
  auto c_d_w_tri_s1_gen_f = [&](const Seqd& s, unsigned int parameter){ return general_eval::cputime_ms_of_method(s, parameter, capla_tri_s1_f); };
  LineGenerator c_d_w_tri_s1_gen = { c_d_w_tri_s1_gen_f, "Tri s1 double window pla" };
  plot::plot_lines_generated(dataset
      , {15, 30, 45, 60, 75, 90}
      , {
	paa_gen
	, apca_gen
	, pla_gen
	, d_w_apla_gen
	, d_w_proj_apla_gen
	, c_d_w_mean_gen
	, c_d_w_mean_s1_gen
	, c_d_w_tri_gen
	, c_d_w_tri_s1_gen
      }, p);
  }



  

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


  //RandomWalk walk( NormalFunctor(1) ); 
  //walk.gen_steps(30);
  //walk.save_walk("./tsv/testwalk1.tsv");
  vector<double> dataset2 = parse_tsv("./tsv/testwalk1.tsv",-1);
  Series s = { dataset2, "Normal Walk" };
  // plot_series(s, "./img/");

  vector<double> test = { 7, 5, 5, 3, 2, 4 };
  //auto apca_test = apca::apca(test, 3);

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
