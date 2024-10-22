#include <string>
#include <vector>
#include "./tsv_parsing.h"


namespace ucr_parsing {

  enum DatasetType { TEST, TRAIN, TEST_APPEND_TRAIN };

  std::vector<float> parse_ucr_dataset(std::string dataset_name, std::string dataset_loc, DatasetType type);

  std::vector<std::string> parse_folder_names(std::string directory_path);
};
