#ifndef SEQUENTIAL_SCAN_H
#define SEQUENTIAL_SCAN_H

#include <vector>
#include <iostream>


namespace seq_scan {
double l2_sqr(const double* const s1start, const double* const s2start, unsigned int len);

std::vector<unsigned int> find_similar_subseq_indexes(const std::vector<double>& series, const std::vector<double>& query, double epsilon);
std::vector<unsigned int> find_k_closest_indexes(const std::vector<double>& series, const std::vector<double>& query, unsigned int k);
};

#endif
