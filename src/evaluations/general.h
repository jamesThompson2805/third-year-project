#ifndef EVAL_GENERAL_H
#define EVAL_GENERAL_H

#include "pla.h"
#include <functional>
#include <vector>

/**
 * @file general.h
 * @brief Header file declaring many functions for generating insights into approximations
 */

/**
 * @brief type definition of a DRT, a function that takes a series and number of parameters and returns the uncompressed approximation that used that number of parameters.
 */
typedef std::function<Seqd(const Seqd&, unsigned int)> DRT;
/**
 * @brief type definition of a comparison function, taking two series and returning some distance measure between the two (examples maybe L2 or maximum deviation).
 */
typedef std::function<double(const Seqd&, const Seqd&)> COMP;

/**
 * @brief general_eval is a namespace containing generic methods to compare a series with approximations of itself.
 */
namespace general_eval {
  /**
   * @brief comp_of_method compares a series with the uncompressed approximation of itself under the provided distance measure.
   * @param s is a reference to the series used.
   * @param num_params is the number of parameters you wish to grant the DRT to approximate with.
   * @param f is the DRT to be used (eg. PAA).
   * @param comp is the comparison function (eg. L2 or even non-norm such as Squared Error).
   * @return The observed distance between the exact series and its approximation.
   */
  double comp_of_method(const Seqd& s, unsigned int num_params, DRT f, COMP comp);
  /**
   * @brief mse_of_method compares a series with the uncompressed approximation of itself under the mean squared error measure.
   * @param s is a reference to the series used.
   * @param num_params is the number of parameters you wish to grant the DRT to approximate with.
   * @param f is the DRT to be used (eg. PAA).
   * @return The observed distance between the exact series and its approximation.
   */
  double mse_of_method(const Seqd& s, unsigned int num_params, DRT f);
  /**
   * @brief l2_of_method compares a series with the uncompressed approximation of itself under the euclidean distance.
   * @param s is a reference to the series used.
   * @param num_params is the number of parameters you wish to grant the DRT to approximate with.
   * @param f is the DRT to be used (eg. PAA).
   * @return The observed distance between the exact series and its approximation.
   */
  double l2_of_method(const Seqd& s, unsigned int num_params, DRT f);
  /**
   * @brief maxdev_of_method compares a series with the uncompressed approximation of itself under the maximum deviation distance.
   * @param s is a reference to the series used.
   * @param num_params is the number of parameters you wish to grant the DRT to approximate with.
   * @param f is the DRT to be used (eg. PAA).
   * @return The observed distance between the exact series and its approximation.
   */
  double maxdev_of_method(const Seqd& s, unsigned int num_params, DRT f);
  /**
   * @brief cputime_of_method returns the time taken to compute the approximation of the series.
   * @param s is a reference to the series used.
   * @param num_params is the number of parameters you wish to grant the DRT to approximate with.
   * @param f is the DRT to be used (eg. PAA).
   * @return The time taken to form the approximation (the function times an application of f).
   */
  double cputime_ms_of_method(const Seqd& s, unsigned int num_params, DRT f);

  /**
   * @brief get_comp_over_sizes compares a series with the uncompressed approximation of itself for a variety of sizes under the provided distance measure.
   * @param s is a reference to the series used.
   * @param sizes is a vector of sizes the dataset is to be shortened to. Order doesn't matter as copies of the sequence are made but sizes must be smaller than the original.
   * @param num_params is the number of parameters you wish to grant the DRT to approximate with.
   * @param f is the DRT to be used (eg. PAA).
   * @param comp is the comparison function (eg. L2 or even non-norm such as Squared Error).
   * @return The observed distances between the exact series (at various lengths) and its approximation of those lengths.
   */
  Seqd get_comp_over_sizes(const Seqd& s, Sequi sizes, unsigned int num_params, DRT f, COMP comp);
  /**
   * @brief get_mse_over_sizes compares a series with the uncompressed approximation of itself for a variety of sizes under the Mean Squared Error.
   * @param s is a reference to the series used.
   * @param sizes is a vector of sizes the dataset is to be shortened to. Order doesn't matter as copies of the sequence are made but sizes must be smaller than the original.
   * @param num_params is the number of parameters you wish to grant the DRT to approximate with.
   * @param f is the DRT to be used (eg. PAA).
   * @return The observed distances between the exact series (at various lengths) and its approximation of those lengths.
   */
  Seqd get_mse_over_sizes(const Seqd& s, Sequi sizes, unsigned int num_params, DRT f);
  /**
   * @brief get_mse_over_sizes compares a series with the uncompressed approximation of itself for a variety of sizes under the Euclidean distance.
   * @param s is a reference to the series used.
   * @param sizes is a vector of sizes the dataset is to be shortened to. Order doesn't matter as copies of the sequence are made but sizes must be smaller than the original.
   * @param num_params is the number of parameters you wish to grant the DRT to approximate with.
   * @param f is the DRT to be used (eg. PAA).
   * @return The observed distances between the exact series (at various lengths) and its approximation of those lengths.
   */
  Seqd get_l2_over_sizes(const Seqd& s, Sequi sizes, unsigned int num_params, DRT f);
  /**
   * @brief get_comp_over_sizes compares a series with the uncompressed approximation of itself for a variety of sizes under maximum deviation distance.
   * @param s is a reference to the series used.
   * @param sizes is a vector of sizes the dataset is to be shortened to. Order doesn't matter as copies of the sequence are made but sizes must be smaller than the original.
   * @param num_params is the number of parameters you wish to grant the DRT to approximate with.
   * @param f is the DRT to be used (eg. PAA).
   * @return The observed distances between the exact series (at various lengths) and its approximation of those lengths.
   */
  Seqd get_maxdev_over_sizes(const Seqd& s, Sequi sizes, unsigned int num_params, DRT f);
  /**
   * @brief get_comp_over_sizes compares the time to execute the DRT approximation under various sizes
   * @param s is a reference to the series used.
   * @param sizes is a vector of sizes the dataset is to be shortened to. Order doesn't matter as copies of the sequence are made but sizes must be smaller than the original.
   * @param num_params is the number of parameters you wish to grant the DRT to approximate with.
   * @param f is the DRT to be used (eg. PAA).
   * @return The observed distances between the exact series (at various lengths) and its approximation of those lengths.
   */
  Seqd get_cputime_over_sizes(const Seqd& s, Sequi sizes, unsigned int num_params, DRT f);

  /**
   * @brief get_comp_over_num_params compares a series with the uncompressed approximation of itself for a variety of sizes under the provided distance measure.
   * @param s is a reference to the series used.
   * @param sizes is a vector of sizes the dataset is to be shortened to. Order doesn't matter as copies of the sequence are made but sizes must be smaller than the original.
   * @param num_params is the number of parameters you wish to grant the DRT to approximate with.
   * @param f is the DRT to be used (eg. PAA).
   * @param comp is the comparison function (eg. L2 or even non-norm such as Squared Error).
   * @return The observed distances between the exact series (at various lengths) and its approximation of those lengths.
   */
  Seqd get_comp_over_num_params(const Seqd& s, Sequi vec_params, DRT f, COMP comp);
  Seqd get_mse_over_num_params(const Seqd& s, Sequi vec_params, DRT f);
  Seqd get_l2_over_num_params(const Seqd& s, Sequi vec_params, DRT f);
  Seqd get_maxdev_over_num_params(const Seqd& s, Sequi vec_params, DRT f);
  Seqd get_cputime_over_num_params(const Seqd& s, Sequi vec_params, DRT f);

  Seqd get_comp_over_DRTs(const Seqd& s, unsigned int num_params, std::vector<DRT> fs, COMP comp);
  Seqd get_mse_over_DRTs(const Seqd& s, unsigned int num_params, std::vector<DRT> fs);
  Seqd get_l2_over_DRTs(const Seqd& s, unsigned int num_params, std::vector<DRT> fs);
  Seqd get_maxdev_over_DRTs(const Seqd& s, unsigned int num_params, std::vector<DRT> fs);
  Seqd get_cputime_over_DRTs(const Seqd& s, unsigned int num_params, std::vector<DRT> fs);



}


#endif
