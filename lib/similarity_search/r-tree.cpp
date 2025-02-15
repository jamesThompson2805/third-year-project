#include "r-tree.h"

#include <queue>
using std::queue;

template <typename R>
R& RTree<R>::get_mbr(RTreeNodePtr<R> r)
{
  if (r.index() == 0) return std::get<0>(r)->mbr;
  return std::get<1>(r)->mbr;
}
template <typename R>
void RTree<R>::set_mbr(RTreeNodePtr<R> r, R mbr)
{
  if (r.index() == 0) {
    std::get<0>(r)->mbr = mbr;
  } else {
    std::get<1>(r)->mbr = mbr;
  }
}
template <typename R>
void RTree<R>::append_ptr(RTreeNodePtr<R> n, RTreeNodePtr<R> ptr)
{
  if (n.index() == 0) {
    std::get<0>(n)->ptrs.push_back(ptr);
  } else {
    std::get<1>(n)->ptrs.push_back(ptr);
  }
}
template <typename R>
void RTree<R>::remove_ptr_i(RTreeNodePtr<R> n, int ptr_i)
{
  if (n.index() == 0) {
    std::get<0>(n)->ptrs.erase(ptr_i);
  } else {
    std::get<1>(n)->ptrs.erase(ptr_i);
  }
}

template <typename R>
RTree<R>::RTree(unsigned int max_entries, unsigned int min_entries, FPtrArea<R> area_f, FPtrAreaMerge<R> merge_f, FPtrMBRDist<R> p_dist_f)
  : max_entries(max_entries), min_entries(min_entries), area_f(area_f), merge_f(merge_f), uncompr_point_dist_f(p_dist_f), root(nullptr)
{}

template <typename R>
RTree<R>::~RTree()
{
  destroy_tree_from(root);
  root = nullptr;
}

// non recursive destruction to avoid crashing the stack
template <typename R>
void RTree<R>::destroy_tree_from(RTreeNodePtr<R> node)
{
  if (node.index() == 1) { // leaf node
    auto leaf = std::get<1>(node);
    if (leaf == nullptr) return;
    delete leaf;
  } else {
    auto nonleaf = std::get<0>(node);
    if (nonleaf == nullptr) return;
    queue<RTreeNodePtr<R>> q = {node};

    while (q.size() > 0) {
      RTreeNodePtr<R> next = q.pop();
      if (next.index() == 0) { // next is a leaf
	destroy_tree_from(next);
	continue;
      }
      auto nonleaf = std::get<0>(next);
      for (RTreeNodePtr<R> ptr : nonleaf->ptrs) {
	q.push(ptr);
      }
      delete next;
    }
  }
}
