#include "general.h"

#include "pla.h"
#include "mse.h"

#include <cmath>
#include <algorithm>
#include <chrono>
#include <ratio>


double general_eval::comp_of_method(const Seqd& s, unsigned int num_params, DRT f, COMP comp)
{
  Seqd reduced = f(s, num_params);
  return comp(s, reduced);
}
double general_eval::mse_of_method(const Seqd &s, unsigned int num_params, DRT f)
{
  return general_eval::comp_of_method(s, num_params, f, mse::mse_between_seq);
}
double general_eval::maxdev_of_method(const Seqd &s, unsigned int num_params, DRT f)
{
  return general_eval::comp_of_method(s, num_params, f, mse::maxdev_between_seq);
}
double general_eval::cputime_ms_of_method(const Seqd &s, unsigned int num_params, DRT f)
{ // https://www.learncpp.com/cpp-tutorial/timing-your-code/ tutorial for how to time cpp code
  auto start = std::chrono::high_resolution_clock::now();
  Seqd reduced = f(s, num_params);
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  return duration.count();
}


Seqd general_eval::get_comp_over_sizes(const Seqd& s, Sequi sizes, unsigned int num_params, DRT f, COMP comp)
{
  Seqd comparisons;
  Seqd resized;
  for (const auto& ui : sizes) {
    resized.resize( ui );
    std::copy_n(s.begin(), ui, resized.begin());
    Seqd reduced = f(resized, num_params);
    comparisons.push_back( comp(resized, reduced) );
  }
  return comparisons;
}
Seqd general_eval::get_mse_over_sizes(const Seqd& s, Sequi sizes, unsigned int num_params, DRT f)
{
  return general_eval::get_comp_over_sizes(s, sizes, num_params, f, mse::mse_between_seq);
}
Seqd general_eval::get_maxdev_over_sizes(const Seqd& s, Sequi sizes, unsigned int num_params, DRT f)
{
  return general_eval::get_comp_over_sizes(s, sizes, num_params, f, mse::maxdev_between_seq);
}
Seqd get_cputime_over_sizes(const Seqd& s, Sequi sizes, unsigned int num_params, DRT f)
{
  Seqd timings;
  Seqd resized;
  for (const auto& ui : sizes) {
    resized.resize( ui );
    std::copy_n(s.begin(), ui, resized.begin());
    timings.push_back( general_eval::cputime_ms_of_method(resized, num_params, f) );
  }
  return timings;
}


Seqd general_eval::get_comp_over_num_params(const Seqd& s, Sequi vec_params, DRT f, COMP comp)
{
  Seqd comparisons;
  for (const auto& ui : vec_params) {
    Seqd reduced = f(s, ui);
    comparisons.push_back( comp(s, reduced) );
  }
  return comparisons;
}
Seqd general_eval::get_mse_over_num_params(const Seqd& s, Sequi vec_params, DRT f)
{
  return general_eval::get_comp_over_num_params(s, vec_params, f, mse::mse_between_seq);
}
Seqd general_eval::get_maxdev_over_num_params(const Seqd& s, Sequi vec_params, DRT f)
{
  return general_eval::get_comp_over_num_params(s, vec_params, f, mse::maxdev_between_seq);
}
Seqd get_cputime_over_num_params(const Seqd& s, Sequi vec_params, DRT f) {

  Seqd timings;
  for (const auto& ui : vec_params) {
    timings.push_back( general_eval::cputime_ms_of_method(s, ui, f) );
  }
  return timings;
}
