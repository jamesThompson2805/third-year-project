#include "capla.h"

#include <cmath>
#include <algorithm>
#include <functional>
using std::vector;

#include "pla.h"
#include "mse.h"
#include "conv_double_window.h"

double capla_eval::mse_of_method(const Seq& s, unsigned int num_params, const Seq& ldist, const Seq& rdist)
{
  Seq capla = pla::apla_to_seq( c_d_w::conv_pla(s, num_params, ldist, rdist) );
  return std::sqrt( mse::se_between_seq(s, capla) );
}
double capla_eval::maxdev_of_method(const Seq& s, unsigned int num_params, const Seq& ldist, const Seq& rdist)
{
  Seq capla = pla::apla_to_seq( c_d_w::conv_pla(s, num_params, ldist, rdist) );
  return mse::maxdev_between_seq(s, capla);
}


Seq get_comp_over_sizes(const Seq& s, std::vector<unsigned int> sizes, unsigned int num_params, const Seq& ldist, const Seq& rdist, std::function<double(const Seq&, const Seq&)> comparison)
{
  Seq comparisons;
  Seq resized;
  for (const auto& ui : sizes) {
    resized.resize( ui );
    std::copy_n(s.begin(), ui, resized.begin());
    Seq capla = pla::apla_to_seq(c_d_w::conv_pla(resized, num_params, ldist, rdist));
    comparisons.push_back( comparison(resized, capla) );
  }
  return comparisons;
}
Seq capla_eval::get_mse_over_sizes(const Seq& s, std::vector<unsigned int> sizes, unsigned int num_params, const Seq& ldist, const Seq& rdist)
{
  Seq a = get_comp_over_sizes(s, sizes, num_params, ldist, rdist, mse::se_between_seq);
  std::for_each(a.begin(), a.end(), [](auto& d){ d=std::sqrt(d); });
  return a;
}
Seq capla_eval::get_maxdev_over_sizes(const Seq& s, std::vector<unsigned int> sizes, unsigned int num_params, const Seq& ldist, const Seq& rdist)
{
  return get_comp_over_sizes(s, sizes, num_params, ldist, rdist, mse::maxdev_between_seq);
}


Seq get_comp_over_num_params(const Seq& s, std::vector<unsigned int> vec_params, const Seq& ldist, const Seq& rdist, std::function<double(const Seq&, const Seq&)> comp)
{
  Seq comparisons;
  for (const auto& ui : vec_params) {
    Seq capla = pla::apla_to_seq(c_d_w::conv_pla(s, ui, ldist, rdist));
    comparisons.push_back( comp(s, capla) );
  }
  return comparisons;
}
Seq get_mse_over_num_params(const Seq& s, std::vector<unsigned int> vec_params, const Seq& ldist, const Seq& rdist)
{
  Seq a = get_comp_over_num_params(s, vec_params, ldist, rdist, mse::se_between_seq);
  std::for_each(a.begin(), a.end(), [](auto& d){ d=std::sqrt(d); });
  return a;
}
Seq get_maxdev_over_num_params(const Seq& s, std::vector<unsigned int> vec_params, const Seq& ldist, const Seq& rdist)
{
  return get_comp_over_num_params(s, vec_params, ldist, rdist, mse::maxdev_between_seq);
}
