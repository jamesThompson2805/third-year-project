#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "ucr_parsing.h"
#include "paa.h"
#include "pla.h"
#include "z_norm.h"

#include "plotting/series_plotting.h"

#include "plotting/plot_dimreduct_paa.cpp"

using std::vector, std::string;

int main()
{
  string ucr_datasets_loc = "external/data/UCRArchive_2018/";
  vector<string> datasets = ucr_parsing::parse_folder_names("external/data/UCRArchive_2018/");

  std::cout << datasets[0] << std::endl;

  vector<float> x;
  for (int j=1; j<=30; ++j) x.push_back(j);

  vector<Line> lines;
  for (int i=0; i<20; ++i) {
    vector<float> dataset_i_mse;
    vector<float> dataset_i = ucr_parsing::parse_ucr_dataset(datasets[i], ucr_datasets_loc, ucr_parsing::DatasetType::TEST_APPEND_TRAIN);
    z_norm::z_normalise(dataset_i);
    for (int j=1; j<=30; ++j) {
      dataset_i_mse.push_back( paa::paa_mse(dataset_i, j) );
    }
    lines.push_back( { x, dataset_i_mse, "dataset " + std::to_string(i) } );
  }
  PlotLabels pl = {"comparison of paa on first 20 datasets", "interval size", "MSE error"};
  plot_lines( lines, pl);

  return 0;
}
