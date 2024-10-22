#include <fstream>
#include <boost/algorithm/string.hpp>

#include "tsv_parsing.h"

using namespace tsv_parsing;
using boost::algorithm::split;

vector<float> parse_tsv(std::string filename, int max_lines=-1)
{
  std::ifstream ifs(filename);
  std::string line;
  vector<std::string> str_floats;
  int curr_line = 0;

  vector<float> parsed;

  while (getline(ifs, line)) {

    if ( (curr_line >= max_lines) && max_lines > -1) break;
    
    split(str_floats, line, boost::is_any_of("\t"));

    for (auto str : str_floats) {
      if (str == "NaN") continue; // skip Nan values
      parsed.emplace_back( std::stof(str) );
    }
    curr_line++;
  }

  return parsed;
}
