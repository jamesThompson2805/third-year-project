#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "ucr_parsing.h"
#include "z_norm.h"
#include "mse.h"

#include "paa.h"
#include "pla.h"
#include "double_window.h"
#include "dac_curve_fitting.h"
#include "exact_dp.h"

#include "plotting/series_plotting.h"
#include "plotting/plot_dimreduct_paa.cpp"
#include "plotting/plot_dimreduct_pla.cpp"

#include "random_walk.h"

using std::vector, std::string;

int main()
{
  using namespace ucr_parsing;

  string ucr_datasets_loc = "external/data/UCRArchive_2018/";
  vector<string> datasets = parse_folder_names(ucr_datasets_loc);

  /*
  for (int i=0; i<datasets.size(); i++) {
    std::cout << i << ": " << datasets[i] << "\t";
    if ((i+1)%3 == 0) std::cout << std::endl;
  }
  */

  unsigned int di = 25;
  vector<double> dataset = parse_ucr_dataset(datasets[di], ucr_datasets_loc,  DatasetType::TRAIN);
  dataset.resize(240);
  z_norm::z_normalise(dataset);


  auto paa_f = [](const vector<double>& s, unsigned int num_params){ return paa::paa_to_seq(paa::paa(s, num_params), (s.size() / num_params) + (s.size() % num_params != 0)); };
  auto pla_f = [](const vector<double>& s, unsigned int num_params){ return pla::pla_to_seq(pla::pla(s, num_params), (2*s.size() / num_params) + (2*s.size() % num_params != 0)); };

  auto d_w_apca_f = [](const vector<double>& s, unsigned int num_params){ return paa::apca_to_seq( d_w::simple_paa(s, num_params, 5, 5)); };
  auto d_w_proj_apca_f = [](const vector<double>& s, unsigned int num_params){ return paa::apca_to_seq( d_w::y_proj_paa(s, num_params, 5, 5)); };

  auto d_w_apla_f = [](const vector<double>& s, unsigned int num_params){ return pla::apla_to_seq( d_w::simple_pla(s, num_params, 5, 5)); };
  auto d_w_proj_apla_f = [](const vector<double>& s, unsigned int num_params){ return pla::apla_to_seq( d_w::y_proj_pla(s, num_params, 5, 5)); };

  auto exact_apaa_f = [](const vector<double>& s, unsigned int num_params){ return paa::apca_to_seq( exact_dp::exact_paa(s, num_params) ); }; 
  auto exact_apla_f = [](const vector<double>& s, unsigned int num_params){ return pla::apla_to_seq( exact_dp::exact_pla(s, num_params) ); }; 

  //auto apaa = exact_dp::exact_paa(dataset, 30);
  // plot_any_apaa_subseq( dataset, datasets[di], [](const vector<double>& s){ return exact_dp::exact_paa(s, 30); }, 0, 119);
  plot_any_apla_subseq( dataset, datasets[di], [](const vector<double>& s){ return exact_dp::exact_pla(s, 30); }, 0, 119);
  
  return 0;
}
