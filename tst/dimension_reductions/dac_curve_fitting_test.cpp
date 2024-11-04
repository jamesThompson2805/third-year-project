#include "dac_curve_fitting.h"
#include <gtest/gtest.h>
#include <vector>

TEST(DAC_TEST, DAC_CONSTANT_SINGLETON) {
  std::vector<double> f = { 3 };
  std::vector<std::tuple<double, unsigned int>> res = { { 3.0, 0 } };
  EXPECT_EQ( dac_curve_fitting::dac_constant(f, 1), res);
}
  
TEST(DAC_TEST, DAC_CONSTANT_EXACT) {
  std::vector<double> f = { 3, 3, 3, 4, 4 };
  std::vector<std::tuple<double, unsigned int>> res = { { 3.0, 2 }, {4.0, 4} };
  EXPECT_EQ( dac_curve_fitting::dac_constant(f, 0), res);
}

TEST(DAC_TEST, DAC_CONSTANT_SMALL_DEVIATIONS) {
  std::vector<double> f = { 2.5, 3, 3.5, 10.5, 9.5 };
  std::vector<std::tuple<double, unsigned int>> res = { { 3.0, 2 }, {10.0, 4} };
  EXPECT_EQ( dac_curve_fitting::dac_constant(f, 0.5), res);
}
