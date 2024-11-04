#ifndef UCR_PARSING_H
#define UCR_PARSING_H

#include <string>
#include <vector>

namespace ucr_parsing {

  std::vector<double> parse_tsv(std::string filename, int max_lines);

  enum DatasetType { TEST, TRAIN, TEST_APPEND_TRAIN };

  std::vector<double> parse_ucr_dataset(std::string dataset_name, std::string dataset_loc, DatasetType type);

  std::vector<std::string> parse_folder_names(std::string directory_path);
};

#endif
