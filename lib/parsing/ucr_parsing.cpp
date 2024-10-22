#include "ucr_parsing.h"
using namespace ucr_parsing;

#include <filesystem>
#include <algorithm>

std::vector<float> ucr_parsing::parse_ucr_dataset(std::string dataset_name, std::string dataset_loc, DatasetType type)
{
  std::string dataset_type;
  if (type == DatasetType::TEST) {
    return tsv_parsing::parse_tsv(dataset_loc + dataset_name + "/" + dataset_name + "_TEST.tsv");
  } else if (type == DatasetType::TRAIN) {
    return tsv_parsing::parse_tsv( dataset_loc + dataset_name + "/" + dataset_name + "_TRAIN.tsv");
  } else {
    std::vector<float> test = tsv_parsing::parse_tsv( dataset_loc + dataset_name + "/" + dataset_name + "_TEST.tsv");
    std::vector<float> train = tsv_parsing::parse_tsv( dataset_loc + dataset_name + "/" + dataset_name + "_TRAIN.tsv");
    test.insert( test.end(), train.begin(), train.end());
    return train;
  }
}

std::vector<std::string> ucr_parsing::parse_folder_names(std::string directory_path)
{
  std::vector<std::string> datasets;
  for (const auto& entry : std::filesystem::directory_iterator( directory_path )) {
    datasets.emplace_back(entry.path().filename());
  }
  std::sort(datasets.begin(), datasets.end());
  return datasets;
}
