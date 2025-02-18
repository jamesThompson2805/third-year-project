#ifndef R_TREE_H
#define R_TREE_H

#include <vector>
#include <array>
#include <set>
#include <variant>

#include <queue>
using std::queue;

template <typename R, typename I>
struct RTreeNode;

template <typename R, typename I>
struct LeafEntry {
  R mbr;
  I st_index;
  LeafEntry(R mbr, I st_index) : mbr(mbr), st_index(st_index) {}
};
template <typename R, typename I>
struct RTreeNode {
  R mbr;
  RTreeNode<R,I>* parent;
  std::variant<std::vector<LeafEntry<R,I>>, std::vector<RTreeNode<R,I>*>> entries;
  RTreeNode(R mbr, RTreeNode<R,I>* p, std::variant<std::vector<LeafEntry<R,I>>, std::vector<RTreeNode<R,I>*>> e) : mbr(mbr), parent(p), entries(e) {}
};

template <typename R>
using FPtrArea = double (*)(const R&);
template <typename R>
using FPtrAreaMerge = R (*)(const R&, const R&);
template <typename R>
using FPtrMBRDistSqr = double (*)(const std::vector<double>& q, const R&);


#include <iostream>
template <typename R, typename I>
class RTree {
private:
  RTreeNode<R,I>* root;
  FPtrArea<R> area_f;
  FPtrAreaMerge<R> merge_f;
  FPtrMBRDistSqr<R> uncompr_point_dist_sqr_f;
  unsigned int max_entries;
  unsigned int min_entries;

public:
  RTree(unsigned int max_entries, unsigned int min_entries, FPtrArea<R> area_f, FPtrAreaMerge<R> merge_f, FPtrMBRDistSqr<R> p_dist_f)
    : max_entries(max_entries), min_entries(min_entries), area_f(area_f), merge_f(merge_f), uncompr_point_dist_sqr_f(p_dist_f), root(nullptr) {}
  ~RTree()
  {
    destroy_tree_from(root);
    root = nullptr;
  }

  void insert(R mbr, unsigned int st_index);

  std::vector<I> sim_search(const std::vector<double>& q, double epsilon);


  RTreeNode<R,I>* split_node(RTreeNode<R,I>* const node);

  std::array<unsigned int, 2> quad_pick_seeds(const std::vector<const R*>& mbrs);
  std::array<unsigned int, 2> lin_pick_seeds(const RTreeNode<R,I>* const node);
  unsigned int quad_pick_next_node(const std::vector<const R*>& mbrs, const std::set<unsigned int>& g1, const R& g1mbr, const std::set<unsigned int>& g2, const R& g2mbr);
  unsigned int lin_pick_next_node(const std::vector<const R*>& mbrs, const std::set<unsigned int>& g1, const R& g1mbr, const std::set<unsigned int>& g2, const R& g2mbr);

private:

  RTreeNode<R,I>* choose_leaf(const R&);
  RTreeNode<R,I>* choose_leaf_from(const R&, RTreeNode<R,I>*);

  void adjust_tree(RTreeNode<R,I>* node); 

  std::vector<const R*> get_entry_mbrs( const RTreeNode<R,I> *const n);
  R rebuild_mbr( const RTreeNode<R,I> *const n);
  unsigned int get_num_mbrs( const RTreeNode<R,I> *const n);

  void destroy_tree_from(RTreeNode<R,I>* node);
};

/* ******************** General R Tree functions ************************ */
// non recursive destruction to avoid crashing the stack
template <typename R, typename I>
void RTree<R,I>::destroy_tree_from(RTreeNode<R,I>* node)
{
  if (node == nullptr) return;
  if (node->entries.index() == 0) { // leaf node
    delete node;
  } else {
    queue<RTreeNode<R,I>*> q;
    q.push(node);
    std::set<RTreeNode<R,I>*> seen;

    while (q.size() > 0) {
      RTreeNode<R,I>* next = q.front();
      q.pop();
      seen.insert(next);
      if (next->entries.index() == 0) { // next is a leaf
	destroy_tree_from(next);
	continue;
      }
      for (RTreeNode<R,I>* entry : std::get<1>(next->entries)) {
	q.push(entry);
      }
      delete next;
    }
  }
  stopNow: return;
}


template <typename R, typename I>
std::vector<const R*> RTree<R,I>::get_entry_mbrs( const RTreeNode<R,I> *const n)
{
  if (n->entries.index() == 0) {
    std::vector<const R*> mbrs;
    for ( const LeafEntry<R,I>& l : std::get<0>(n->entries) ){
      mbrs.push_back(&l.mbr);
    }
    return mbrs;
  } else {
    std::vector<const R*> mbrs;
    for ( const RTreeNode<R,I>* ep : std::get<1>(n->entries) ){
      mbrs.push_back(&ep->mbr);
    }
    return mbrs;
  }
}

template <typename R, typename I>
unsigned int RTree<R,I>::get_num_mbrs( const RTreeNode<R,I> *const n)
{
  if (n->entries.index() == 0)
    return std::get<0>(n->entries).size();
  return std::get<1>(n->entries).size();
}

template <typename R, typename I>
R RTree<R,I>::rebuild_mbr( const RTreeNode<R,I> *const n)
{
  std::vector<const R*> mbrs = get_entry_mbrs(n);
  R mbr = *mbrs[0];
  for (const R* const rp : mbrs) {
    mbr = merge_f( mbr, *rp );
  }
  return mbr;
}

/* *********************************** Insert definitions ******************* */
#include <cmath>

template <typename R, typename I>
RTreeNode<R,I>* RTree<R,I>::choose_leaf(const R& mbr)
{
  return choose_leaf_from(mbr, root);
}

// recursive call here is all good as goes directly downwards, doesn't branch at all
template <typename R, typename I>
RTreeNode<R,I>* RTree<R,I>::choose_leaf_from(const R& mbr, RTreeNode<R,I>* node)
{
  // if root is leaf node
  if (node->entries.index() == 0) return node;
  
  // otherwise find the child rect that increases area least to accomodate it and recursively call
  double min_area_incr = -1;
  RTreeNode<R,I>* min_r;
  for (RTreeNode<R,I>* e_ptr : std::get<1>(node->entries)) {
    double area_incr = area_f( merge_f( mbr, e_ptr->mbr ) ) - area_f(e_ptr->mbr);
    if (min_area_incr < 0 || area_incr < min_area_incr || area_incr == min_area_incr && area_f(e_ptr->mbr) < area_f(min_r->mbr)) {
      min_r = e_ptr;
      min_area_incr = area_incr;
    }
  }
  return choose_leaf_from(mbr, min_r);
}

template <typename R, typename I>
std::array<unsigned int, 2> RTree<R,I>::quad_pick_seeds(const std::vector<const R*>& mbrs)
{
  double max_d = -1;
  unsigned int i1, i2;

  for (unsigned int i=0; i<mbrs.size(); i++) {
    for (unsigned int j=0; j<i; j++) {
      double d = area_f( merge_f(*mbrs[i], *mbrs[j]) )
		  - area_f( *mbrs[i] )
		  - area_f( *mbrs[j] );
      if (d > max_d) {
	max_d = d;
	i1 = i;
	i2 = j;
      }
    }
  }
  return { i1, i2 };
}

template <typename R, typename I>
unsigned int RTree<R,I>::quad_pick_next_node(const std::vector<const R*>& mbrs, const std::set<unsigned int>& g1, const R& g1mbr, const std::set<unsigned int>& g2, const R& g2mbr)
{
  double max_d = 0;
  int max_i = -1;

  for (int i=0; i<mbrs.size(); i++) {
    if ( g1.count(i) > 0 || g2.count(i) > 0) continue;
    double d = area_f( merge_f(g1mbr, *mbrs[i]) )
		- area_f( merge_f(g2mbr, *mbrs[i]) );
    if (std::abs(d) >= std::abs(max_d)) {
      max_d = d;
      max_i = i;
    }
  }
  return max_i;
}

// split node only adjusts the MBR of node and new twin_node
//  does not adjust mbr of parent
//  may invalidate max entries of parent
//  must be called with adjust tree to maintain invariant
template <typename R, typename I>
RTreeNode<R,I>* RTree<R,I>::split_node(RTreeNode<R,I>* node)
{
  std::vector<const R*> mbrs = get_entry_mbrs(node);

  std::array<unsigned int, 2> seeds = quad_pick_seeds(mbrs);
  std::set<unsigned int> g0 = {seeds[0]}, g1 = {seeds[1]};
  R g0mbr = *mbrs[seeds[0]], g1mbr = *mbrs[seeds[1]];

  int num_remaining = get_num_mbrs(node) - 2;

  while (g0.size() < max_entries && g1.size() < max_entries && g0.size()+num_remaining > min_entries && g1.size()+num_remaining > min_entries) {
    unsigned int next_i = quad_pick_next_node(mbrs, g0, g0mbr, g1, g1mbr);
    if (next_i == -1) break; // exhausted all numbers
    num_remaining--;
    double area_incr_g0 = area_f( merge_f( g0mbr, *mbrs[next_i] ) ) - area_f( g0mbr );
    double area_incr_g1 = area_f( merge_f( g1mbr, *mbrs[next_i] ) ) - area_f( g1mbr );
    if (area_incr_g0 < area_incr_g1 || area_incr_g0 == area_incr_g1 && area_f(g0mbr) < area_f(g1mbr)) {
      g0.insert(next_i);
      g0mbr = merge_f(g0mbr, *mbrs[next_i]);
    } else {
      g1.insert(next_i);
      g1mbr = merge_f(g1mbr, *mbrs[next_i]);
    }
  }
  if (g0.size()==max_entries || g1.size()+num_remaining <= min_entries) {
    unsigned int next_i = quad_pick_next_node(mbrs, g0, g0mbr, g1, g1mbr);
    while (next_i != -1) {
      g1.insert(next_i);
      g1mbr = merge_f(g1mbr, *mbrs[next_i]);
      next_i = quad_pick_next_node(mbrs, g0, g0mbr, g1, g1mbr);
    }
  }
  if (g1.size()==max_entries || g0.size()+num_remaining <= min_entries) {
    unsigned int next_i = quad_pick_next_node(mbrs, g0, g0mbr, g1, g1mbr);
    while (next_i != -1) {
      g0.insert(next_i);
      g0mbr = merge_f(g0mbr, *mbrs[next_i]);
      next_i = quad_pick_next_node(mbrs, g0, g0mbr, g1, g1mbr);
    }
  }

  if ( node->entries.index() == 0 ) {
    std::vector<LeafEntry<R,I>>& leaves = std::get<0>(node->entries);

    std::vector<LeafEntry<R,I>> g0entries;
    std::vector<LeafEntry<R,I>> g1entries;
    for (unsigned int n : g1) {
      g1entries.push_back( leaves[n] );
    }
    for (unsigned int n : g0) {
      g0entries.push_back( leaves[n] );
    }

    // cover is root case
    if (node->parent == nullptr) {
      RTreeNode<R,I>* new_root = new RTreeNode<R,I>({ node->mbr, nullptr, std::vector<RTreeNode<R,I>*>{node} });
      node->parent = new_root;
      root = new_root;
    }

    node->entries = g0entries;
    node->mbr = g0mbr; // adjust node's mbr to reflect its current entries
    RTreeNode<R,I>* g1node = new RTreeNode<R,I>({g1mbr, node->parent, g1entries});
    std::get<1>(node->parent->entries).push_back( g1node ); // add g1node to parent

    return g1node;

  } else {
    const std::vector<RTreeNode<R,I>*>& entries = std::get<1>(node->entries);

    std::vector<RTreeNode<R,I>*> g0entries;
    std::vector<RTreeNode<R,I>*> g1entries;
    for (unsigned int n : g1) {
      g1entries.push_back(entries[n]);
    } // copy across pointers
    for (unsigned int n : g0) {
      g0entries.push_back(entries[n]);
    } // copy across pointers

    // cover is root case
    if (node->parent == nullptr) {
      RTreeNode<R,I>* new_root = new RTreeNode<R,I>({ node->mbr, nullptr, std::vector<RTreeNode<R,I>*>{node} });
      node->parent = new_root;
      root = new_root;
    }

    RTreeNode<R,I>* g1node = new RTreeNode<R,I>({g1mbr, node->parent, g1entries});
    for (unsigned int n : g1) {
	entries[n]->parent = g1node;
    } // change child entries to have new parent
    node->entries = g0entries;
    node->mbr = g0mbr; // adjust node's mbr to reflect its current entriesk

    std::get<1>(node->parent->entries).push_back( g1node ); // add g1node to parent
    return g1node;
  }
}

template <typename R, typename I>
void RTree<R,I>::adjust_tree(RTreeNode<R,I>* node)
{
  RTreeNode<R,I>* curr_node = node;
  while (curr_node != root) {
    RTreeNode<R,I>* parent_node = curr_node->parent;
    parent_node->mbr = rebuild_mbr(parent_node);
    if (get_num_mbrs(parent_node) > max_entries) {
      split_node(parent_node);
    }
    curr_node = parent_node;
  }
}

template <typename R, typename I>
void RTree<R,I>::insert(R mbr, unsigned int st_index)
{
  if (root == nullptr) {
    root = new RTreeNode<R,I>(mbr, nullptr, std::vector<LeafEntry<R,I>>({{mbr,st_index}}) );
    return;
  }
  if (st_index % 1000 == 0)
    std::cout << "st: " << st_index << " " << std::endl;
  RTreeNode<R,I>* target_leaf = choose_leaf(mbr);
  std::vector<LeafEntry<R,I>>& entries = std::get<0>(target_leaf->entries);
  entries.push_back({mbr, st_index});
  target_leaf->mbr = merge_f( target_leaf->mbr, mbr );

  RTreeNode<R,I>* twin_node = nullptr;
  if (entries.size() > max_entries) {
    twin_node = split_node(target_leaf);
  }
  adjust_tree(target_leaf);
    
}
/* *************************** R Tree Search methods ************************* */

template <typename R, typename I>
std::vector<I> RTree<R,I>::sim_search(const std::vector<double>& q, double epsilon)
{
  epsilon = epsilon * epsilon;
  if (uncompr_point_dist_sqr_f(q, root->mbr) > epsilon) {
    std::cout << root->mbr[0] << " " << root->mbr[1] << std::endl;
    for ( RTreeNode<R,I>* n : std::get<1>(root->entries) ) {
      std::cout <<"	" << n->mbr[0] << " " << n->mbr[1] << std::endl;
    }

    return {};
  }

  queue<const RTreeNode<R,I>*> nodes;
  nodes.push(root);
  std::vector<I> results;

  while (nodes.size() != 0) {
    const RTreeNode<R,I>* const next = nodes.front();
    nodes.pop();
    if (next->entries.index() == 0) {
      for (const LeafEntry<R,I>& l : std::get<0>(next->entries)) {
	if (uncompr_point_dist_sqr_f(q, l.mbr) <= epsilon)
	  results.push_back(l.st_index);
      }
    } else {
      for (const RTreeNode<R,I>* const n : std::get<1>(next->entries)) {
	if (uncompr_point_dist_sqr_f(q, n->mbr) <= epsilon)
	  nodes.push(n);
      }
    }
  }
  return results;
}


#endif
