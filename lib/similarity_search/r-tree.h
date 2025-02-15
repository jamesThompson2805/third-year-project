#ifndef R_TREE_H
#define R_TREE_H

#include <vector>
#include <array>
#include <tuple>
#include <variant>

template <typename R>
struct RTreeNode;
template <typename R>
struct RTreeLeafNode;

template <typename R>
using RTreeNodePtr = std::variant<RTreeNode<R>*, RTreeLeafNode<R>*>;

template <typename R>
struct LeafEntry {
  R mbr;
  unsigned int st_index;
};
template <typename R>
struct RTreeLeafNode {
  R mbr;
  RTreeNodePtr<R> parent;
  std::vector<LeafEntry<R>> st_ptrs;
};
template <typename R>
struct RTreeNode {
  R mbr;
  RTreeNodePtr<R> parent;
  std::vector<RTreeNodePtr<R>> ptrs;
};

template <typename R>
using FPtrArea = double (*)(const R&);
template <typename R>
using FPtrAreaMerge = R (*)(const R&, const R&);
template <typename R>
using FPtrMBRDist = double (*)(const std::vector<double>& q, const R&);


template <typename R>
class RTree {
private:
  RTreeNodePtr<R> root;
  FPtrArea<R> area_f;
  FPtrAreaMerge<R> merge_f;
  FPtrMBRDist<R> uncompr_point_dist_f;
  unsigned int max_entries;
  unsigned int min_entries;

public:
  RTree(unsigned int max_entries, unsigned int min_entries, FPtrArea<R> area_f, FPtrAreaMerge<R> merge_f, FPtrMBRDist<R> p_dist_f);
  ~RTree();

  void insert(R mbr, unsigned int st_index);
  std::vector<unsigned int> search(const R& q_mbr);

  std::vector<unsigned int> sim_search(const std::vector<double>& q, double epsilon);
  std::vector<unsigned int> knn_search(const std::vector<double>& q);


  RTreeNodePtr<R> split_node(RTreeNodePtr<R> node);

  std::array<RTreeNodePtr<R>, 2> quad_pick_seeds_leaf(RTreeLeafNode<R>* node);
  std::array<RTreeNodePtr<R>, 2> quad_pick_seeds_nonleaf(RTreeNode<R>* node);
  std::array<std::vector<RTreeNodePtr<R>>, 2> pick_seeds_lin(RTreeNodePtr<R> node);
  void quad_pick_next_node_leaf(RTreeLeafNode<R>* node, std::array<RTreeLeafNode<R>*, 2>& groups);
  void quad_pick_next_node_nonleaf(RTreeNodePtr<R> node, std::array<RTreeNodePtr<R>, 2>& groups);

private:

  // helper functions for std variant
  R& get_mbr(RTreeNodePtr<R> r);
  void set_mbr(RTreeNodePtr<R> r, R mbr);
  void append_ptr(RTreeNodePtr<R> n, RTreeNodePtr<R> ptr);
  void remove_ptr_i(RTreeNodePtr<R> n, int ptr_i);

  RTreeLeafNode<R>* choose_leaf(const R&);
  RTreeLeafNode<R>* choose_leaf_from(const R&, const RTreeNodePtr<R>&);

  void destroy_tree_from(RTreeNodePtr<R> node);



};


#endif
