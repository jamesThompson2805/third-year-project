#ifndef APLA_SEGMENT_AND_MERGE_H
#define APLA_SEGMENT_AND_MERGE_H

#include "pla.h"
#include <vector>
#include <tuple>

typedef std::vector<std::tuple<DoublePair, unsigned int>> Seqddt;
typedef std::vector<double> Seqd;

namespace segmerge {
  void merge_1(const Seqd& s, Seqddt& s_compr);
  void merge_k(const Seqd& s, Seqddt& s_compr, unsigned int k);
  void merge_to_dim(const Seqd& s, Seqddt& s_compr, unsigned int k);
  void segment_1(const Seqd& s, Seqddt& s_compr);
  void segment_k(const Seqd& s, Seqddt& s_compr, unsigned int k);
  void segment_to_dim(const Seqd& s, Seqddt& s_compr, unsigned int k);
  
  void maxdev_merge_1(const Seqd& s, Seqddt& s_compr);
  void maxdev_merge_k(const Seqd& s, Seqddt& s_compr, unsigned int k);
  void maxdev_segment_1(const Seqd& s, Seqddt& s_compr);
  void maxdev_segment_k(const Seqd& s, Seqddt& s_compr, unsigned int k);

}

#endif
