#include <vector>
#include <iostream>


namespace seq_scan {
float mse_l2(const float* const s1start, const float* const s2start, unsigned int len);

std::vector<int> find_similar_subseq_indexes(const std::vector<float>& series, const std::vector<float>& query, float epsilon);
};
