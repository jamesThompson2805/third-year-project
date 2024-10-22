#include <iostream>
#include <string>
#include <vector>
using std::vector;

#include "ucr_parsing.h"

int main()
{
  vector<std::string> datasets = ucr_parsing::parse_folder_names("external/data/UCRArchive_2018/");

  std::cout << datasets[0] << std::endl;
  return 0;
}
