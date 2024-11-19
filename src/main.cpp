#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "double_window.h"
#include "ucr_parsing.h"
#include "paa.h"
#include "pla.h"
#include "dac_curve_fitting.h"
#include "z_norm.h"
#include "mse.h"

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
  // dataset.resize(280);
  z_norm::z_normalise(dataset);

  // plot_paa_mse(datasets, ucr_datasets_loc, 0, 10, 100);
  // plot_pla_mse(datasets, ucr_datasets_loc, 0, 10, 100);

  // plot_paa_subseq_series_comparison(dataset, datasets[di], 20, 0, 100);
  // plot_pla_subseq_series_comparison(dataset, datasets[di], 20, 0, 100);
  // plot_apaa_subseq_series_comparison(dataset, datasets[di], 0.00005, 0, 100);
  // plot_apla_subseq_series_comparison(dataset, datasets[100], 0.00005, 0, 100);

  // /*
  dataset.resize(120);
  plot_paa_subseq(dataset, datasets[di], 30, 0, 119);
  plot_pla_subseq(dataset, datasets[di], 30, 0, 119);
  std::cout << "PAA MSE: " 
	    << mse::mse_between_seq(dataset,
		paa::paa_to_seq(paa::paa(dataset, 30),4))
	    << std::endl;
  std::cout << "PLA MSE: " 
	    << mse::mse_between_seq(dataset,
		pla::pla_to_seq(pla::chunk_regression(dataset, 30),8))
	    << std::endl;
  // */

  // /*
  plot_any_apaa_subseq( dataset, datasets[di], [](const vector<double>& s){ return d_w::y_proj_paa(s, 30, 2,2); }, 0, 119);
  plot_any_apla_subseq( dataset, datasets[di], [](const vector<double>& s){ return d_w::y_proj_pla(s, 30, 2,2); }, 0, 119);
  std::cout << "APCA MSE: " 
	    << mse::mse_between_seq(dataset,
		paa::apca_to_seq(d_w::y_proj_paa(dataset, 30, 2, 2)))
	    << std::endl;
  std::cout << "APLA MSE: " 
	    << mse::mse_between_seq(dataset,
		pla::apla_to_seq(d_w::y_proj_pla(dataset, 30, 2, 2)))
	    << std::endl;

  // */

  // /*
  plot_any_apaa_subseq( dataset, datasets[di] + " for mean of window", [](const vector<double>& s){ return d_w::simple_paa(s, 30, 2,2); }, 0, 119);
  plot_any_apla_subseq( dataset, datasets[di] + " for mean of window", [](const vector<double>& s){ return d_w::simple_pla(s, 30, 2,2); }, 0, 119);
  std::cout << "APCA MSE: " 
	    << mse::mse_between_seq(dataset,
		paa::apca_to_seq(d_w::simple_paa(dataset, 30, 2, 2)))
	    << std::endl;
  std::cout << "APLA MSE: " 
	    << mse::mse_between_seq(dataset,
		pla::apla_to_seq(d_w::simple_pla(dataset, 30, 2, 2)))
	    << std::endl;
  // */


  return 0;
}


