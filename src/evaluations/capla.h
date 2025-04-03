#ifndef EVAL_CAPLA_H
#define EVAL_CAPLA_H
/**
 * @file capla.h
 * @brief Header file for helpful functions for utilising the ChangeInSlope sliding window technique.
 */


#include "general.h"

#include "pla.h"

/**
 * @brief type definition for a function that takes a series and number of parameters and returns a valid Adaptive PLA approximation.
 */
using DRT_COMPR = std::function<std::vector<std::tuple<DoublePair, unsigned int>> (const Seqd&,unsigned int)>;

/**
 * @brief capla_eval is the namespace holding functions that return functions given parameters for how to construct the 'windows'.
 * For more details about CAPLA, please read the paper also found in the GitHub repository.
 */
namespace capla_eval {
  /**
   * Function creates non-weighted distribution of given window size and returns function to apply to data.
   * @param win_size is the desired size of the windows, 0 is invalid here
   * @return Function that accepts series and number of parameters to return approximation already uncompressed as a sequence of the same size.
   */
  DRT generate_mean_DRT(unsigned int win_size);
  /**
   * Function creates non-weighted distribution of given window size and returns function to apply to data, this time returning in compressed format.
   * @param win_size is the desired size of the windows, 0 is invalid here
   * @return Function that accepts series and number of parameters to return approximation
   */
  DRT_COMPR generate_mean_DRT_COMPR(unsigned int win_size);
  /**
   * Function creates non-weighted distribution of given window size and returns function to apply to data. This is formed of the additional hypothesis that gradient between the two potential lines shouldn't be accounted for.
   * @param win_size is the desired size of the windows, 0 is invalid here
   * @return Function that accepts series and number of parameters to return approximation already uncompressed as a sequence of the same size.
   */
  DRT generate_mean_skip_one_DRT(unsigned int win_size);
  /**
   * Function creates weighted distribution of given window size in a triangular fashion, those closer to the split have greater weight and returns function to apply to data.
   * @param win_size is the desired size of the windows, 0 is invalid here
   * @return Function that accepts series and number of parameters to return approximation already uncompressed as a sequence of the same size.
   */
  DRT generate_tri_DRT(unsigned int win_size);
  /**
   * Function creates weighted distribution of given window size in a triangular fashion, those closer to the split have greater weight and returns function to apply to data. This also has the additional hypothesis that gradient between the two potential lines shouldn't be accounted for.
   * @param win_size is the desired size of the windows, 0 is invalid here
   * @return Function that accepts series and number of parameters to return approximation already uncompressed as a sequence of the same size.
   */
  DRT generate_tri_skip_one_DRT(unsigned int win_size);
};


#endif
