#include "r_tree.h"

#include <gtest/gtest.h>

#include <cmath>
#include <array>
#include <vector>
#include <set>
#include <algorithm>

typedef std::vector<double> MBR1D;

double area_1d(const MBR1D& mbr)
{
  return mbr[1] - mbr[0];
}

MBR1D merge_1d(const MBR1D& a, const MBR1D& b)
{
  return { std::min(a[0], b[0]), std::max(a[1], b[1]) };
}

double mbr_1d_point_dist(const std::vector<double>& q, const MBR1D& r)
{
  if (q.size() == 0) return 0.0;
  if (q[0] < r[0]) return (r[0] - q[0])*(r[0] - q[0]);
  if (q[0] > r[1]) return (r[1] - q[0])*(r[1] - q[0]);
  return 0.0;
}


TEST(RTree, RTreeInsert) {
  std::vector<double> nums;
  RTree<MBR1D, unsigned int> rtree(40, 10, area_1d, merge_1d, mbr_1d_point_dist); 
  
  for (double i=0.0; i<=100'000; i+=0.5) {
    nums.push_back(i);
  }

  for (int i=0; i<nums.size(); i++) {
    rtree.insert({nums[i], nums[i]}, i);
  }

  std::cout << "finished inserting" << std::endl;

  std::vector<double> q = { 10'000 };
  std::vector<unsigned int> searched = rtree.sim_search(q, 1.0);

  for ( unsigned int i : searched ) {
    std::cout << "index " << i << " result " << nums[i] << " ";
  }
  std::cout << std::endl;

  std::set<double> expected_els = { 10'000, 9999.5, 10'000.5, 9999, 10001 };
  std::set<double> result;
  std::for_each( searched.begin(), searched.end(), [&result, &nums](unsigned int i){ result.insert( nums[i] );});

  auto retrieve = [](const unsigned int& n, const std::vector<double>& s){ return std::vector<std::array<const double*,2>>({{&s[n], &s[n]}}); };


  for (auto& [lptr,rptr] : rtree.knn_search({10'000}, 10, retrieve, nums) ) {
    std::cout << *lptr << " ";
  }
  std::cout << std::endl;
  std::cout << rtree.pruning_power({10'000.3}, retrieve, nums) << std::endl;


  EXPECT_EQ(expected_els, result);
}

TEST(RTree, RTreeAPCA) {
  std::vector<double> nums;
  RTree<MBR1D, unsigned int> rtree(40, 10, area_1d, merge_1d, mbr_1d_point_dist); 
  
  for (double i=0.0; i<=100'000; i+=0.5) {
    nums.push_back(i);
  }

  for (int i=0; i<nums.size(); i++) {
    rtree.insert({nums[i], nums[i]}, i);
  }

  std::cout << "finished inserting" << std::endl;

  std::vector<double> q = { 10'000 };
  std::vector<unsigned int> searched = rtree.sim_search(q, 1.0);

  for ( unsigned int i : searched ) {
    std::cout << "index " << i << " result " << nums[i] << " ";
  }
  std::cout << std::endl;

  std::set<double> expected_els = { 10'000, 9999.5, 10'000.5, 9999, 10001 };
  std::set<double> result;
  std::for_each( searched.begin(), searched.end(), [&result, &nums](unsigned int i){ result.insert( nums[i] );});

  auto retrieve = [](const unsigned int& n, const std::vector<double>& s){ return std::vector<std::array<const double*,2>>({{&s[n], &s[n]}}); };


  for (auto& [lptr,rptr] : rtree.knn_search({10'000}, 10, retrieve, nums) ) {
    std::cout << *lptr << " ";
  }
  std::cout << std::endl;
  std::cout << rtree.pruning_power({10'000.3}, retrieve, nums) << std::endl;


  EXPECT_EQ(expected_els, result);
}
