#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <tuple>
#include <algorithm>

#include "ucr_parsing.h"
#include "z_norm.h"
#include "mse.h"

#include "paa.h"
#include "apca.h"

#include "pla.h"
#include "double_window.h"
#include "dac_curve_fitting.h"
#include "exact_dp.h"

#include "conv_double_window.h"

#include "plotting/series_plotting.h"

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

  /*
  unsigned int di = 25;
  vector<double> dataset = parse_ucr_dataset(datasets[di], ucr_datasets_loc,  DatasetType::TRAIN);
  dataset.resize(480);
  z_norm::z_normalise(dataset);
  */


  auto paa_f = [](const vector<double>& s, unsigned int num_params){ return paa::paa(s, num_params), (s.size() / num_params) + (s.size() % num_params != 0); };
  auto pla_f = [](const vector<double>& s, unsigned int num_params){ return pla::pla(s, num_params), (2*s.size() / num_params) + (2*s.size() % num_params != 0); };

  auto d_w_apca_f = [](const vector<double>& s, unsigned int num_params){ return d_w::simple_paa(s, num_params, 5, 5); };
  auto d_w_proj_apca_f = [](const vector<double>& s, unsigned int num_params){ return d_w::y_proj_paa(s, num_params, 5, 5); };

  auto d_w_apla_f = [](const vector<double>& s, unsigned int num_params){ return d_w::simple_pla(s, num_params, 3, 3); };
  auto d_w_proj_apla_f = [](const vector<double>& s, unsigned int num_params){ return d_w::y_proj_pla(s, num_params, 5, 5); };

  auto exact_apaa_f = [](const vector<double>& s, unsigned int num_params){ return exact_dp::min_mse_paa(s, num_params); }; 
  auto exact_apla_f = [](const vector<double>& s, unsigned int num_params){ return exact_dp::min_mse_pla(s, num_params); }; 

  vector<double> ldist = { 1.0/3.0, 1.0/3.0, 1.0/3.0};
  vector<double> rdist = { 0.0, 1.0/2.0, 1.0/2.0};
  auto conv_apla_f = [&ldist, &rdist](const vector<double>& s, unsigned int num_params){ return c_d_w::conv_pla(s, num_params, ldist, rdist); }; 

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

  vector<double> test = { 7, 5, 5, 3, 2, 4, 4, 6 };
  vector<tuple<double, unsigned int>> apca = apca::apca(test, 6);
  for (auto [v,d] : apca) {
    std::cout << v << " : " << d << std::endl;
  }

  //RandomWalk walk( NormalFunctor(1) ); 
  //walk.gen_steps(30);
  //walk.save_walk("./tsv/testwalk1.tsv");
  vector<double> dataset2 = parse_tsv("./tsv/testwalk1.tsv",-1);
  Series s = { dataset2, "Normal Walk" };
  // plot_series(s, "./img/");

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
