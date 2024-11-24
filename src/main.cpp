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
  dataset.resize(240);
  z_norm::z_normalise(dataset);


  // plot_paa_mse(datasets, ucr_datasets_loc, 0, 10, 100);
  // plot_pla_mse(datasets, ucr_datasets_loc, 0, 10, 100);

  // plot_paa_subseq_series_comparison(dataset, datasets[di], 20, 0, 100);
  // plot_pla_subseq_series_comparison(dataset, datasets[di], 20, 0, 100);
  // plot_apaa_subseq_series_comparison(dataset, datasets[di], 0.00005, 0, 100);
  // plot_apla_subseq_series_comparison(dataset, datasets[100], 0.00005, 0, 100);

  auto paa_f = [](const vector<double>& s, unsigned int num_params){ return paa::paa_to_seq(paa::paa(s, num_params), (s.size() / num_params) + (s.size() % num_params != 0)); };
  auto pla_f = [](const vector<double>& s, unsigned int num_params){ return pla::pla_to_seq(pla::pla(s, num_params), (2*s.size() / num_params) + (2*s.size() % num_params != 0)); };

  auto d_w_apca_f = [](const vector<double>& s, unsigned int num_params){ return paa::apca_to_seq( d_w::simple_paa(s, num_params, 5, 5)); };
  auto d_w_proj_apca_f = [](const vector<double>& s, unsigned int num_params){ return paa::apca_to_seq( d_w::y_proj_paa(s, num_params, 5, 5)); };

  auto d_w_apla_f = [](const vector<double>& s, unsigned int num_params){ return pla::apla_to_seq( d_w::simple_pla(s, num_params, 5, 5)); };
  auto d_w_proj_apla_f = [](const vector<double>& s, unsigned int num_params){ return pla::apla_to_seq( d_w::y_proj_pla(s, num_params, 5, 5)); };

  for (int i=3; i<=10; i++) {
    vector<double> x;
    for (int j=4; j<=4000; j+=100) x.push_back(j);

    vector<Line> lines;
    vector<vector<double>> datasets_mse(6);
    vector<double> dataset_i = parse_ucr_dataset(datasets[i], ucr_datasets_loc,  DatasetType::TRAIN);
    for (int j=4; j<=4000; j+=100) {
      datasets_mse[0].push_back( mse::maxdev_between_seq(dataset_i, paa_f(dataset_i, j)) );
      datasets_mse[1].push_back( mse::maxdev_between_seq(dataset_i, pla_f(dataset_i, j)) );
      datasets_mse[2].push_back( mse::maxdev_between_seq(dataset_i, d_w_apca_f(dataset_i, j)) );
      datasets_mse[3].push_back( mse::maxdev_between_seq(dataset_i, d_w_proj_apca_f(dataset_i, j)) );
      datasets_mse[4].push_back( mse::maxdev_between_seq(dataset_i, d_w_apla_f(dataset_i, j)) );
      datasets_mse[5].push_back( mse::maxdev_between_seq(dataset_i, d_w_proj_apla_f(dataset_i, j)) );
    }
    lines.push_back( { x, datasets_mse[0], "paa" } );
    lines.push_back( { x, datasets_mse[1], "pla" } );
    lines.push_back( { x, datasets_mse[2], "sliding window paa" } );
    lines.push_back( { x, datasets_mse[3], "sliding window intervals paa" } );
    lines.push_back( { x, datasets_mse[4], "sliding window pla" } );
    lines.push_back( { x, datasets_mse[5], "sliding window intervals pla" } );
    PlotLabels pl = {"comparison of Maximum Deviation for DRT on dataset" + std::to_string(i), "number of parameters", "MaxDev error"};
    plot_lines( lines, pl);
  }

  /*
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
  */

  /*
  plot_any_apaa_subseq( dataset, datasets[di], [](const vector<double>& s){ return d_w::y_proj_paa(s, 30, 2,2); }, 0, 119);
  plot_any_apla_subseq( dataset, datasets[di], [](const vector<double>& s){ return d_w::y_proj_pla(s, 30, 2,2); }, 0, 119);
  std::cout << "APCA Y PROJECTION MSE: " 
	    << mse::mse_between_seq(dataset,
		paa::apca_to_seq(d_w::y_proj_paa(dataset, 30, 2, 2)))
	    << std::endl;
  std::cout << "APLA Y PROJECTION MSE: " 
	    << mse::mse_between_seq(dataset,
		pla::apla_to_seq(d_w::y_proj_pla(dataset, 30, 2, 2)))
	    << std::endl;

  */

  /*
  plot_any_apaa_subseq( dataset, datasets[di] + " for mean of window", [](const vector<double>& s){ return d_w::simple_paa(s, 30, 2,2); }, 0, 119);
  plot_any_apla_subseq( dataset, datasets[di] + " for mean of window", [](const vector<double>& s){ return d_w::simple_pla(s, 30, 2,2); }, 0, 119);
  std::cout << "APCA MEAN MSE: " 
	    << mse::mse_between_seq(dataset,
		paa::apca_to_seq(d_w::simple_paa(dataset, 30, 2, 2)))
	    << std::endl;
  std::cout << "APLA MEAN MSE: " 
	    << mse::mse_between_seq(dataset,
		pla::apla_to_seq(d_w::simple_pla(dataset, 30, 2, 2)))
	    << std::endl;
  */


  return 0;
}
