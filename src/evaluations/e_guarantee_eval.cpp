#include "e_guarantee_eval.h"

#include <chrono>



double precision_eval::compr_ratio_of_method(const Seqd& s, double epsilon, DRT_COMPR_PG f)
{
  return (double)f(s,epsilon).size() / (double)s.size();
}
double precision_eval::cputime_of_method(const Seqd& s, double epsilon, DRT_COMPR_PG f)
{
  auto start = std::chrono::high_resolution_clock::now();
  Seqddt approx = f(s, epsilon);
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  return duration.count();
}

Seqd precision_eval::get_segments_over_epsilons(const Seqd& s, Seqd epsilons, DRT_COMPR_PG f)
{
  Seqd ratios;
  for (auto e : epsilons) {
    ratios.push_back( precision_eval::compr_ratio_of_method(s, e, f) );
  }
  return ratios;
}

Seqd precision_eval::get_cputime_over_epsilons(const Seqd& s, Seqd epsilons, DRT_COMPR_PG f)
{
  Seqd timings;
  for (auto e : epsilons) {
    timings.push_back( precision_eval::cputime_of_method(s, e, f) );
  }
  return timings;
}
