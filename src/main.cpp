#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

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

  vector<double> dataset = parse_ucr_dataset(datasets[100], ucr_datasets_loc,  DatasetType::TEST);
  z_norm::z_normalise(dataset);

  plot_paa_mse(datasets, ucr_datasets_loc, 0, 10, 100);
  plot_pla_mse(datasets, ucr_datasets_loc, 0, 10, 100);

  plot_paa_subseq_series_comparison(dataset, datasets[100], 3, 0, 100);
  plot_pla_subseq_series_comparison(dataset, datasets[100], 3, 0, 100);
  plot_apaa_subseq_series_comparison(dataset, datasets[100], 0.00005, 0, 100);
  plot_apla_subseq_series_comparison(dataset, datasets[100], 0.00005, 0, 100);

  
  return 0;
}


