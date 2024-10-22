#include <string>
#include <vector>
using std::vector;


namespace tsv_parsing {
  vector<float> parse_tsv(std::string filename, int max_lines=-1);
};
