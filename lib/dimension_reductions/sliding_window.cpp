#include "sliding_window.h"

Seqddt sw::sliding_window(const Seqd &q, double epsilon)
{
  Seqddt segments;
  unsigned int start_i = 0;
  unsigned int end_i = 0;
  while (end_i < q.size()) {
    DoublePair seg_approx = pla::regression(q.data() + start_i, q.data() + end_i);
    bool within_epsilon = true;
    for (int i=start_i; i<=end_i; i++) {
      within_epsilon &= (std::abs( q[i] - seg_approx[0] - seg_approx[1]*(i-start_i) ) < epsilon);
    }
    if (within_epsilon) {
      end_i++;
      continue;
    }
    segments.push_back({ pla::regression(q.data()+start_i, q.data()+end_i-1), start_i });
    start_i = end_i;
  }
  return segments;
}
