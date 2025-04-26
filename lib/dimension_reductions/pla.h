#ifndef PLA_H
#define PLA_H

#include <array>
#include <vector>
#include <tuple>

#include <functional>
#include <algorithm>

/**
 * @file pla.h contains the pla approximation and utilitiy definitions adjacent to it
 */

/**
 * @brief DoublePair describes a line using two parameters, the first is the y-intercept and the second the gradient
 * It is customary for a PLA algorithm to have each segments linear approximation approximate as if the segment starts at 0
 */
typedef std::array<double, 2> DoublePair;

/**
 * @brief Seqd is an array of double values, representing a sequence
 */
typedef std::vector<double> Seqd;
/**
 * @brief Seqdd is an array of pair of doubles, representing a series of lines
 */
typedef std::vector<DoublePair> Seqdd;
/**
 * @brief Sequi is an array of unsigned integers
 */
typedef std::vector<unsigned int> Sequi;
/**
 * @brief Seqddt is an array of triples, representing the information for an adaptive PLA segment
 */
typedef std::vector<std::tuple<DoublePair,unsigned int>> Seqddt;

/**
 * @brief pla namespace holds functions for PLA algorithm and decompression
 */
namespace pla {
  /**
   * @brief chunk_series takes a series and returns the series segmented into chunks
   * @param series is the series to segment into equal sized pieces
   * @param chunk_size is the size of the pieces
   * @return array of sequences representing the pieces
   */
std::vector<std::vector<double>> chunk_series(std::vector<double> series, unsigned int chunk_size);


/**
 * @brief regression calculates the line of best fit on an array of points
 * @param start points to first element
 * @param end points to last element
 * @return line of best fit of elements
 */
DoublePair regression(const double* const start, const double* const end);
/**
 * @brief regression_thru_point calculates line of best fit that starts at first point and approximates rest
 * @param first_point is inserted at beginning of series and line starts from this point
 * @param start is first element
 * @param end points to last element
 * @return line of best fit that must go through the first point
 */
double regression_thru_point(double first_point, const double* const start, const double* const end);

/**
 * @brief subsequence_regression calculates the line of best fit for every subsequence of size w
 * @param series is the series to take subsequences from
 * @param w is the size of the subsequences
 * @return an array of all the lines of best fit for every subsequence
 */
std::vector<DoublePair> subsequence_regression( const std::vector<double>& series, unsigned int w);

/**
 * @brief pla calculates the pla algorithm on a series
 * @param series is the series to approximate
 * @param num_params is the target dimension
 * @return a sorted array of lines, each for the corresponding segment
 */
std::vector<DoublePair> pla( const std::vector<double>& series, unsigned int num_params);

/**
 * @brief pla_mse calculates pla on the series and returns the mean squared error of the approximation to the original
 * @param series is the series to approximate and compare to
 * @param interval_size is the size of the interval
 * @return mean squared error of approximation to original
 */
double pla_mse(const std::vector<double>& series, unsigned int interval_size);


/**
 * @brief pla_to_seq takes a pla approximation and its segment size and returns the full size approximation
 * @param pla_s is the compressed series
 * @param int_size is the size of a segment
 * @return a series representing the uncompressed pla_s
 */
std::vector<double> pla_to_seq(const std::vector<DoublePair> pla_s, unsigned int int_size);
/**
 * @brief apla_to_seq takes a apla approximation and returns the uncompressed approximation
 * @param apla_s is the compressed series
 * @return a series representing the uncompressed apla_s
 */
std::vector<double> apla_to_seq(const std::vector< std::tuple<DoublePair, unsigned int>> apla_s);

/**
 * @brief APLA_DRT represents a function that takes a series and target dimension and returns the sorted array of APLA segments
 */
using APLA_DRT = std::function< std::vector<std::tuple<DoublePair, unsigned int>>(const std::vector<double>&, unsigned int)>;
/**
 * @brief APLA is a fixed length array of segments that ensures the dimension of an approximation are all the same
 */
template <unsigned int NS>
using APLA = std::array<std::tuple<DoublePair,unsigned int>,NS>;
/**
 * @brief apla_drt_on_subseqs calculates Adaptive PLA on every subsequence of the given size and returns the array of approximations
 * @param q is the series to have its subsequences approximated
 * @param subseq_size is the size of each subsequence
 * @param f is the Adaptive PLA method to approximate the subsequences
 * @return an array of all the approximations of each subsequence
 */
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
