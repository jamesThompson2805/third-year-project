#ifndef APLA_SEGMENT_AND_MERGE_H
#define APLA_SEGMENT_AND_MERGE_H

#include "pla.h"

/**
 * @file apla_segment_and_merge.h is file containing the split and merge algorithm for Adaptive PLA representations
 */


/**
 * @brief segmerge namespace holds functions for splitting and merging approximations given the original sequence
 */
namespace segmerge {
  /**
   * @brief merge_1 function merges one segment in the approximation with its neighbour, merging only the closest two
   * @param s is the original sequence
   * @param s_compr is its compressed representation
   */
  void merge_1(const Seqd& s, Seqddt& s_compr);
  /**
   * @brief merge_k function merges k segments in the approximation with its neighbour, merging only the closest two each time
   * @param s is the original sequence
   * @param s_compr is its compressed representation
   * @param k is the number of merges to perform
   */
  void merge_k(const Seqd& s, Seqddt& s_compr, unsigned int k);
  /**
   * @brief merge_to_dim merges until the number of segments is k/3, making the dimension of the approximation k
   * @param s is the original sequence
   * @param s_compr is its compressed representation
   * @param k is the dimension to leave the approximation at
   */
  void merge_to_dim(const Seqd& s, Seqddt& s_compr, unsigned int k);
  /**
   * @brief segment_1 function splits one segment in the approximation
   * @param s is the original sequence
   * @param s_compr is its compressed representation
   */
  void segment_1(const Seqd& s, Seqddt& s_compr);
  /**
   * @brief segment_k function splits k segments in the approximation
   * @param s is the original sequence
   * @param s_compr is its compressed representation
   * @param k is the number of segments to perform
   */
  void segment_k(const Seqd& s, Seqddt& s_compr, unsigned int k);
  /**
   * @brief segment_to_dim splits until the number of segments is k/3, making the dimension of the approximation k
   * @param s is the original sequence
   * @param s_compr is its compressed representation
   * @param k is the dimension to leave the approximation at
   */
  void segment_to_dim(const Seqd& s, Seqddt& s_compr, unsigned int k);
  
  void segment_k_opt(const Seqd& s, Seqddt& s_compr, unsigned int k);
  void merge_k_opt(const Seqd& s, Seqddt& s_compr, unsigned int k);
}

#endif
