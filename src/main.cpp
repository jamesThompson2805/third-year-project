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

#include "plotting/series_plotting.h"

#include "plotting/plot_dimreduct_paa.cpp"
#include "plotting/plot_dimreduct_pla.cpp"

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

  vector<double> dataset = parse_ucr_dataset(datasets[44], ucr_datasets_loc,  DatasetType::TRAIN);
  dataset.resize(280);
  z_norm::z_normalise(dataset);
  Series s = {dataset, datasets[44]};
  plot_series(s);

  // plot_paa_mse(datasets, ucr_datasets_loc, 0, 10, 100);
  // plot_pla_mse(datasets, ucr_datasets_loc, 0, 10, 100);

  // plot_paa_subseq_series_comparison(dataset, datasets[100], 20, 0, 100);
  // plot_pla_subseq_series_comparison(dataset, datasets[100], 20, 0, 100);
  // plot_apaa_subseq_series_comparison(dataset, datasets[100], 0.00005, 0, 100);
  // plot_apla_subseq_series_comparison(dataset, datasets[100], 0.00005, 0, 100);

  
  // plot_any_apaa_subseq( dataset, datasets[10], [](const vector<double>& s){ return d_w::simple_paa(s, 30, 2,2); }, 0, 89);
  // plot_any_apla_subseq( dataset, datasets[10], [](const vector<double>& s){ return d_w::y_proj_pla(s, 72, 2,2); }, 0, 119);
  //plot_any_apaa_subseq( dataset, datasets[10]+" for y proj", [](const vector<double>& s){ return d_w::y_proj_paa(s, 80, 10,10); }, 0, 400);
  // plot_paa_subseq(dataset, datasets[10], 20, 0, 100);


  return 0;
}


