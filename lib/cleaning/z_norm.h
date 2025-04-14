#ifndef Z_NORM_H
#define Z_NORM_H

#include <vector>

/**
 * @file z_norm.h
 * @brief Header file containing the means to z-normalise data.
 * The only header file of the data cleaning library as no other methods needed to be implemented in the duration of the project.
 * It offers only the function z_normalise that takes a mutable series reference and mutates all values to have mean zero, variance 1.
 */

/**
 * @brief z_norm namespace contains z_normalisation function.
 */
namespace z_norm {
  /**
   * @brief z_normalise function takes a mutable reference to a sequence and mutates the sequence to be Z Normalised.
   * @param series is a mutable reference to a sequence.
   */
  void z_normalise(std::vector<double>& series);
};

#endif
