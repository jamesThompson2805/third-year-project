#ifndef PLA_H
#define PLA_H

#include <array>
#include <vector>
#include <tuple>

#include <functional>
#include <algorithm>

typedef std::array<double, 2> DoublePair;

typedef std::vector<double> Seqd;
typedef std::vector<DoublePair> Seqdd;
typedef std::vector<unsigned int> Sequi;
typedef std::vector<std::tuple<DoublePair,unsigned int>> Seqddt;

namespace pla {
std::vector<std::vector<double>> chunk_series(std::vector<double> series, unsigned int chunk_size);


DoublePair regression(const double* const start, const double* const end);
double regression_thru_point(double first_point, const double* const start, const double* const end);

// w is num items compressed to a linear function
std::vector<DoublePair> sliding_window_regression( const std::vector<double>& series, unsigned int w);

std::vector<DoublePair> chunk_regression( const std::vector<double>& series, unsigned int num_params);
std::vector<DoublePair> pla( const std::vector<double>& series, unsigned int num_params);

double pla_mse(const std::vector<double>& series, unsigned int interval_size);

std::vector<double> pla_to_seq(const std::vector<DoublePair>, unsigned int int_size);
std::vector<double> apla_to_seq(const std::vector< std::tuple<DoublePair, unsigned int>>);

using APLA_DRT = std::function< std::vector<std::tuple<DoublePair, unsigned int>>(const std::vector<double>&, unsigned int)>;
template <unsigned int NS>
using APLA = std::array<std::tuple<DoublePair,unsigned int>,NS>;
template <unsigned int NS>
std::vector<APLA<NS>> apla_drt_on_subseqs( const std::vector<double>& q, unsigned int subseq_size, APLA_DRT f)
{
  std::vector<APLA<NS>> subseqs_compr;
  for (int i=0; i<q.size() - subseq_size; i++) {
    std::vector<double> q_i( q.cbegin()+i, q.cbegin()+i+subseq_size); 
    auto apla = f(q_i,NS*3);
    APLA<NS> apla_arr;
    std::copy(q_i.begin(), q_i.end(), apla_arr.begin());
    subseqs_compr.push_back( apla_arr );
  }
  return subseqs_compr;
}
};

#endif
