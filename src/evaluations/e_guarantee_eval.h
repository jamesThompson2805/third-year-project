#ifndef E_GUARANTEE_EVAL_H
#define E_GUARANTEE_EVAL_H

#include "pla.h"
#include <functional>
#include <vector>

/**
 * @file e_guarantee_eval.h
 * @brief Header file declaring functions for analysing PLA algorithms that operate via an epsilon precision guarantee
 */ 

/**
 * @brief type definition of a epsilon precision DRT, taking series and epsilon and returns the compressed approximation that used that has the guarantee
 */
typedef std::function<Seqddt(const Seqd&, double)> DRT_COMPR_PG;


/**
 * @brief precision_eval Namespace encompasses functions for evaluating on epsilon guarantee algorithms
 */
namespace precision_eval {
  
  /**
   * @brief compr_ratio_of_method returns the compression ratio of DRT with the precision.
   * @param s is a reference to the series used.
   * @param epsilon is the epsilon for the precision guarantee.
   * @param f is the DRT to be used (eg. Top Down).
   * @return The compression ratio.
   */
  double compr_ratio_of_method(const Seqd& s, double epsilon, DRT_COMPR_PG f);
  /**
   * @brief segments_of_method returns the time DRT took to achieve the precision
   * @param s is a reference to the series used.
   * @param epsilon is the epsilon for the precision guarantee.
   * @param f is the DRT to be used (eg. Top Down).
   * @return The time taken by the DRT.
   */
  double cputime_of_method(const Seqd& s, double epsilon, DRT_COMPR_PG f);

  /**
   * @brief get_segments_over_epsilon returns the compression ratios of the DRT with the epsilons.
   * @brief compr_ratio_of_method returns the compression ratio of DRT with the precision.
   * @param s is a reference to the series used.
   * @param epsilons is a vector of epsilons for the precision guarantee.
   * @param f is the DRT to be used (eg. Top Down).
   * @return The compression ratios.
   */
  Seqd get_segments_over_epsilons(const Seqd& s, Seqd epsilons, DRT_COMPR_PG f);

  /**
   * @brief get_cputime_over_epsilon returns the number of segments the DRT used to achieve the precision for a variety of epsilons.
   * @param s is a reference to the series used.
   * @param epsilons is a vector of epsilons for the precision guarantee
   * @param f is the DRT to be used (eg. Top Down).
   * @return The time taken by the DRT for each epsilon.
   */
  Seqd get_cputime_over_epsilons(const Seqd& s, Seqd epsilons, DRT_COMPR_PG f);

};

#endif
