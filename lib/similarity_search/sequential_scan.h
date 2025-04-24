#ifndef SEQUENTIAL_SCAN_H
#define SEQUENTIAL_SCAN_H

#include <vector>

/**
 * @file sequential_scan.h
 * @brief Header file for sequential scan methods for similarity search and K nearest neighbours
 */

/**
 * @brief seq_scan is namespace for solutions to simsearch and k-nn by sequential scan
 */
namespace seq_scan {

  /**
   * @brief l2_sqr returns the squared error of the two sequences
   * @param s1start points to beginning of s1 array
   * @param s2start points to beginning of s2 array
   * @param len is the length of both arrays
   * @return the squared error between the sequences
   */
  double l2_sqr(const double* const s1start, const double* const s2start, unsigned int len);

  /**
   * @brief find_similar_subseq_indexes finds all subsequences of series that are within epsilon of query
   * @param series is the large time series to search for subsequences in
   * @param query is the query sequence to search for similar sequences to
   * @param epsilon is the maximum allowed l2 error between a query and returned subseqence
   * @return array of integers representing the start index of a subsequence within epsilon
   */
  std::vector<unsigned int> find_similar_subseq_indexes(const std::vector<double>& series, const std::vector<double>& query, double epsilon);
  /**
   * @brief find_k_closest_indexes finds all the k closest subsequences to a query
   * @param series is the large time series to search for subsequences in
   * @param query is the query sequence to search for similar sequences to
   * @param k is the number of subsequences to find
   * @return array of integers representing the start index of the k closest subsequences
   */
  std::vector<unsigned int> find_k_closest_indexes(const std::vector<double>& series, const std::vector<double>& query, unsigned int k);
};

#endif
