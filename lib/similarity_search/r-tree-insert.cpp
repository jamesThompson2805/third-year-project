#include "r-tree.h"
#include <algorithm>

template <typename R>
RTreeLeafNode<R>* RTree<R>::choose_leaf(const R& mbr)
{
  return choose_leaf_from(mbr, root);
}

// recursive call here is all good as goes directly downwards, doesn't branch at all
template <typename R>
RTreeLeafNode<R>* RTree<R>::choose_leaf_from(const R& mbr, const RTreeNodePtr<R>& node)
{
  // if root is leaf node
  if (node.index() == 1) return std::get<1>(node);
  
  // otherwise find the child rect that increases area least to accomodate it and recursively call
  double min_area_incr = -1;
  RTreeNodePtr<R> min_r;
  for (RTreeNodePtr<R> e_ptr : std::get<0>(node)->ptrs) {
    double area_incr = area_f( merge_f( get_mbr(min_r), get_mbr(e_ptr) ) ) - area_f(get_mbr(e_ptr));
    if (min_area_incr < 0 || area_incr < min_area_incr || area_incr == min_area_incr && area_f(get_mbr(e_ptr)) < area_f(get_mbr(min_r))) {
      min_r = e_ptr;
      min_area_incr = area_incr;
    }
  }
  return choose_leaf_from(mbr, min_r);
}

template <typename R>
std::array<RTreeNodePtr<R>, 2> RTree<R>::quad_pick_seeds_leaf(RTreeLeafNode<R>* node)
{
  double max_d = -1;
  int i1, i2;

  for (int i=0; i<node->st_ptrs.size(); i++) {
    for (int j=0; j<i; j++) {
      double d = area_f( merge_f(node->st_ptrs[i].mbr, node->st_ptrs[j].mbr) )
		  - area_f( node->st_ptrs[i].mbr )
		  - area_f( node->st_ptrs[j].mbr );
      if (d > max_d) {
	max_d = d;
	i1 = i;
	i2 = j;
      }
    }
  }

  RTreeLeafNode<R>* node1 = new RTreeLeafNode<R>( {node->st_ptrs[i1].mbr, node, {node->st_ptrs[i1].st_index}} );
  RTreeLeafNode<R>* node2 = new RTreeLeafNode<R>( {node->st_ptrs[i2].mbr, node, {node->st_ptrs[i2].st_index}} );

  node->st_ptrs.erase(i1);
  node->st_ptrs.erase(i2);

  return { node1, node2 };
}
template <typename R>
std::array<RTreeNodePtr<R>, 2> RTree<R>::quad_pick_seeds_nonleaf(RTreeNode<R>* node)
{
  double max_d = -1;
  int i1, i2;

  for (int i=0; i<node->ptrs.size(); i++) {
    for (int j=0; j<i; j++) {
      double d = area_f( merge_f( get_mbr(node->ptrs[i]), get_mbr(node->ptrs[j])) )
		  - area_f( get_mbr(node->ptrs[i]) )
		  - area_f( get_mbr(node->ptrs[j]) );
      if (d > max_d) {
	max_d = d;
	i1 = i;
	i2 = j;
      }
    }
  }

  if (node->ptrs[i1].index() == 0) {
    RTreeNode<R>* node1 = new RTreeNode<R>( { node->ptrs[i1].mbr, node, { node->ptrs[i1] } } );
    RTreeNode<R>* node2 = new RTreeNode<R>( { node->ptrs[i2].mbr, node, { node->ptrs[i2] } } );
    node->ptrs[i1].parent = node1;
    node->ptrs[i2].parent = node2;

    node->ptrs.erase(i1);
    node->ptrs.erase(i2);

    return { node1, node2 };
  } else {
    RTreeLeafNode<R>* node1 = new RTreeNode<R>( { node->ptrs[i1].mbr, node, { node->ptrs[i1] } } );
    RTreeLeafNode<R>* node2 = new RTreeNode<R>( { node->ptrs[i2].mbr, node, { node->ptrs[i2] } } );
    node->ptrs[i1].parent = node1;
    node->ptrs[i2].parent = node2;

    node->ptrs.erase(i1);
    node->ptrs.erase(i2);

    return { node1, node2 };
  }
}

template <typename R>
void RTree<R>::quad_pick_next_node_leaf(RTreeLeafNode<R>* node, std::array<RTreeLeafNode<R>*, 2>& groups)
{
  double max_d = 0;
  int max_i;

  for (int i=0; i<node->st_ptrs.size(); i++) {
    double d = area_f( merge_f(groups[0]->mbr, node->st_ptrs[i].mbr) )
		- area_f( merge_f(groups[1]->mbr, node->st_ptrs[i].mbr) );
    if (std::abs(d) > std::abs(max_d)) {
      max_d = d;
      max_i = i;
    }
  }
  if (max_d < 0) {
    groups[0]->mbr = merge_f( groups[0]->mbr, node->st_ptrs[max_i].mbr );
    groups[0]->st_ptrs.push_back( node->st_ptrs[max_i] );
  } else {
    groups[1]->mbr = merge_f( groups[1]->mbr, node->st_ptrs[max_i].mbr );
    groups[1]->st_ptrs.push_back( node->st_ptrs[max_i] );
  }
  node->st_ptrs.erase(max_i);
}
template <typename R>
void quad_pick_next_node_nonleaf(RTreeNodePtr<R> node, std::array<RTreeNodePtr<R>, 2>& groups)
{
  double max_d = 0;
  int max_i;

  for (int i=0; i<node->st_ptrs.size(); i++) {
    double d = area_f( merge_f(get_mbr(groups[0]), node->ptrs[i].mbr) )
		- area_f( merge_f(get_mbr(groups[1]), node->ptrs[i].mbr) );
    if (std::abs(d) > std::abs(max_d)) {
      max_d = d;
      max_i = i;
    }
  }
  if (max_d < 0 || ) {
    set_mbr( groups[0],  merge_f( groups[0]->mbr, node->st_ptrs[max_i].mbr ) );
    append_ptr( groups[0], node->st_ptrs[max_i] ); 
  } else {
    set_mbr( groups[1], merge_f( groups[1]->mbr, node->st_ptrs[max_i].mbr ) );
    append_ptr( groups[1], node->st_ptrs[max_i] );
  }
  remove_ptr_i(node, max_i);
}

template <typename R>
RTreeNodePtr<R> RTree<R>::split_node(RTreeNodePtr<R> node)
{
   
}

template <typename R>
void RTree<R>::insert(R mbr, unsigned int st_index)
{
  if (root == nullptr) {
    root = std::make_shared<RTreeLeafNode>({mbr, nullptr, {{mbr,st_index}}});
  }
  sp<RTreeLeafNode<R>> target_leaf = choose_leaf(mbr);
  if ( target_leaf->st_ptrs.size() < max_entries) {
    target_leaf->st_ptrs.push_back(st_index);
    target_leaf->mbr = merge_f(target_leaf->mbr, mbr);
  } else {
    
  }
}
