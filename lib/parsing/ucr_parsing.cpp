#include "ucr_parsing.h"

#include <filesystem>
#include <algorithm>
#include <fstream>
#include <boost/algorithm/string.hpp>
using boost::algorithm::split;
using std::vector;

vector<double> ucr_parsing::parse_tsv(std::string filename, int max_lines=-1)
{
  std::ifstream ifs(filename);
  std::string line;
  vector<std::string> str_doubles;
  int curr_line = 0;

  vector<double> parsed;

  while (getline(ifs, line)) {

    if ( (curr_line >= max_lines) && max_lines > -1) break;
    
    split(str_doubles, line, boost::is_any_of("\t"));

    for (auto str : str_doubles) {
      if (str == "NaN") continue; // skip Nan values
      if (str == "") continue; // skip empty values
      parsed.emplace_back( std::stod(str) );
    }
    curr_line++;
  }

  return parsed;
}

vector<double> ucr_parsing::parse_ucr_tsv(std::string filename, int max_lines=-1)
{
  std::ifstream ifs(filename);
  std::string line;
  vector<std::string> str_doubles;
  int curr_line = 0;

  vector<double> parsed;

  bool is_first_in_line = true;
  while (getline(ifs, line)) {

    if ( (curr_line >= max_lines) && max_lines > -1) break;
    
    split(str_doubles, line, boost::is_any_of("\t"));

    is_first_in_line = true;
    for (auto str : str_doubles) {
      if (str == "NaN") continue; // skip Nan values
      if (str == "") continue; // skip empty values
      if (is_first_in_line) {
	is_first_in_line = false;
	continue; // skip first value, it denotes the classes of a time series file
      }
      parsed.emplace_back( std::stod(str) );
    }
    curr_line++;
  }

  return parsed;
}
vector<double> ucr_parsing::parse_ucr_dataset(std::string dataset_name, std::string dataset_loc, DatasetType type)
{
  std::string dataset_type;
  if (type == DatasetType::TEST) {
    return parse_ucr_tsv(""+dataset_loc + dataset_name + "/" + dataset_name + "_TEST.tsv");
  } else if (type == DatasetType::TRAIN) {
    return parse_ucr_tsv(""+dataset_loc + dataset_name + "/" + dataset_name + "_TRAIN.tsv");
  } else {
    vector<double> test = parse_ucr_tsv( ""+dataset_loc + dataset_name + "/" + dataset_name + "_TEST.tsv");
    vector<double> train = parse_ucr_tsv( ""+dataset_loc + dataset_name + "/" + dataset_name + "_TRAIN.tsv");
    train.insert( train.end(), test.begin(), test.end());
    return test;
  }
}

vector<std::string> ucr_parsing::parse_folder_names(std::string directory_path)
{
  vector<std::string> datasets;
  for (const auto& entry : std::filesystem::directory_iterator( directory_path )) {
    datasets.emplace_back(entry.path().filename());
  }
  std::sort(datasets.begin(), datasets.end());
  return datasets;
}
