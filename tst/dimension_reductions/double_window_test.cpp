#include "double_window.h"
#include <gtest/gtest.h>
#include <vector>

TEST(D_W_TEST, D_W_TEST_PAA_SINGLE_SPLIT) {
  std::vector<double> f = { 3.0, 4.0 };
  std::vector<std::tuple<double, unsigned int>> res = { {3.0, 0}, {4.0, 1} };
  EXPECT_EQ( d_w::simple_paa(f, 4, 1, 1), res);
}

TEST(D_W_TEST, D_W_TEST_PAA_THREE_TRIPLES_W1) {
  std::vector<double> f = { 3.0, 3.0, 3.0, 5.0, 5.0, 5.0, 7.0, 7.0, 7.0};
  std::vector<std::tuple<double, unsigned int>> res = { {3.0, 2}, {5.0, 5}, {7.0, 8} };
  EXPECT_EQ( d_w::simple_paa(f, 6, 1, 1), res);
}
TEST(D_W_TEST, D_W_TEST_PAA_THREE_TRIPLES_W2) {
  std::vector<double> f = { 3.0, 3.0, 3.0, 5.0, 5.0, 5.0, 7.0, 7.0, 7.0};
  std::vector<std::tuple<double, unsigned int>> res = { {3.0, 2}, {5.0, 5}, {7.0, 8} };
  EXPECT_EQ( d_w::simple_paa(f, 6, 2, 2), res);
}
TEST(D_W_TEST, D_W_TEST_PAA_THREE_TRIPLES_W3) {
  std::vector<double> f = { 3.0, 3.0, 3.0, 5.0, 5.0, 5.0, 7.0, 7.0, 7.0};
  std::vector<std::tuple<double, unsigned int>> res = { {3.0, 2}, {5.0, 5}, {7.0, 8} };
  EXPECT_EQ( d_w::simple_paa(f, 6, 3, 3), res);
}

TEST(D_W_TEST, D_W_TEST_PAA_THREE_INTERVALS_W2) {
  std::vector<double> f = { 3.0, 3.0, 5.0, 5.0, 5.0, 7.0, 7.0, 7.0, 7.0};
  std::vector<std::tuple<double, unsigned int>> res = { {3.0, 1}, {5.0, 4}, {7.0, 8} };
  EXPECT_EQ( d_w::simple_paa(f, 6, 2, 2), res);
}

TEST(D_W_TEST, D_W_TEST_PLA_SINGLE_SPLIT_W1) {
  std::vector<double> f = { 3.0, 3.0, 4.0, 4.0};
  std::vector<std::tuple<DoublePair, unsigned int>> res = { { {3.0, 0.0}, 1}, { {4.0, 0.0}, 3} };
  EXPECT_EQ( d_w::simple_pla(f, 6, 1, 1), res);
}
TEST(D_W_TEST, D_W_TEST_PLA_SINGLE_SPLIT_W2) {
  std::vector<double> f = { 3.0, 3.0, 3.0, 4.0, 4.0};
  std::vector<std::tuple<DoublePair, unsigned int>> res = { { {3.0, 0.0}, 2}, { {4.0, 0.0}, 4} };
  EXPECT_EQ( d_w::simple_pla(f, 6, 2, 2), res);
}
TEST(D_W_TEST, D_W_TEST_PLA_THREE_TRIPLES_W1) {
  std::vector<double> f = { 3.0, 3.0, 3.0, 5.0, 5.0, 5.0, 7.0, 7.0, 7.0};
  std::vector<std::tuple<double, unsigned int>> res = { {{3.0, 2}, {5.0, 5}, {7.0, 8} };
  EXPECT_EQ( d_w::simple_pla(f, 6, 1, 1), res);
}
TEST(D_W_TEST, D_W_TEST_PLA_THREE_TRIPLES_W2) {
  std::vector<double> f = { 3.0, 3.0, 3.0, 5.0, 5.0, 5.0, 7.0, 7.0, 7.0};
  std::vector<std::tuple<double, unsigned int>> res = { {3.0, 2}, {5.0, 5}, {7.0, 8} };
  EXPECT_EQ( d_w::simple_pla(f, 6, 2, 2), res);
}
TEST(D_W_TEST, D_W_TEST_PLA_THREE_TRIPLES_W3) {
  std::vector<double> f = { 3.0, 3.0, 3.0, 5.0, 5.0, 5.0, 7.0, 7.0, 7.0};
  std::vector<std::tuple<double, unsigned int>> res = { {3.0, 2}, {5.0, 5}, {7.0, 8} };
  EXPECT_EQ( d_w::simple_pla(f, 6, 3, 3), res);
}
