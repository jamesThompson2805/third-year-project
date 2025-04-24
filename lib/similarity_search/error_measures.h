#ifndef MSE_H
#define MSE_H

#include <vector>

/**
 * @file error_measures.h
 * @brief Header file for finding distances between sequences in various manners
 */

/**
 * @brief error_measures namespace handles utility functions for distance measurement
 */
namespace error_measures {
  /**
   * @brief se_between_seq returns the squared error between two sequences
   * @param s1 is a constant reference to the first sequence
   * @param s2 is a constant reference to the second sequence
   */
  double se_between_seq(const std::vector<double>& s1, const std::vector<double>& s2);
  /**
   * @brief mse_between_seq returns the mean squared error between two sequences
   * @param s1 is a constant reference to the first sequence
   * @param s2 is a constant reference to the second sequence
   */
  double mse_between_seq(const std::vector<double>& s1, const std::vector<double>& s2);
  /**
   * @brief se_between_ptrs returns the squared error between two sequences, passed as pointers
   * @param s1_start is a constant pointer to a first sequence start
   * @param s1_end is a constant pointer to a first sequence end
   * @param s2_start is a constant pointer to a second sequence start
   * @param s2_end is a constant pointer to a second sequence end
   */
  double se_between_ptrs(const double* const s1_start, const double* const s1_end, const double* const s2_start, const double* const s2_end);
  /**
   * @brief l2_between_seq returns the euclidean distance between two sequences
   * @param s1 is a constant reference to the first sequence
   * @param s2 is a constant reference to the second sequence
   */
  double l2_between_seq(const std::vector<double>& s1, const std::vector<double>& s2);
  /**
   * @brief maxdev_between_seq returns the maximum deviation between two sequences
   * @param s1 is a constant reference to the first sequence
   * @param s2 is a constant reference to the second sequence
   */
  double maxdev_between_seq(const std::vector<double>& s1, const std::vector<double>& s2);
};
#endif
