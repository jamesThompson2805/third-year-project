#ifndef R_TREE_H
#define R_TREE_H

#include <vector>
#include <array>
#include <set>
#include <variant>

#include <queue>
using std::queue;

#include "error_measures.h"

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
template <typename I>
using FPtrRetrievalMethod = std::vector<std::array<const double*, 2>> (*)(const I&, const std::vector<double>&);


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

  unsigned int total_num_entries;

public:
  RTree(unsigned int max_entries, unsigned int min_entries, FPtrArea<R> area_f, FPtrAreaMerge<R> merge_f, FPtrMBRDistSqr<R> p_dist_f)
    : max_entries(max_entries)
      , min_entries(min_entries)
      , area_f(area_f)
      , merge_f(merge_f)
      , uncompr_point_dist_sqr_f(p_dist_f)
      , root(nullptr)
      , total_num_entries(0) {}
  ~RTree()
  {
    destroy_tree_from(root);
    root = nullptr;
  }

  unsigned int get_size_tree();
  unsigned int get_size_leaves();
  inline unsigned int get_num_leaves() { return total_num_entries; }

  void insert(R mbr, I st_index);

  std::vector<I> sim_search(const std::vector<double>& q, double epsilon);
  std::vector<std::array<const double*, 2>> knn_search(const std::vector<double>& q, unsigned int k, FPtrRetrievalMethod<I> retrieve_f, const std::vector<double>& s);
  std::vector<std::array<const double*, 2>> sim_search_exact(const std::vector<double>& q, double epsilon, FPtrRetrievalMethod<I> retrieve_f, const std::vector<double>& s);
  double pruning_power(const std::vector<double>& q, FPtrRetrievalMethod<I> retrieve_f, const std::vector<double>& s);


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
  unsigned int get_size_tree(const RTreeNode<R,I>* node);
  unsigned int get_size_leaves(const RTreeNode<R,I>* node);

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

    while (q.size() > 0) {
      RTreeNode<R,I>* next = q.front();
      q.pop();
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

template <typename R, typename I>
unsigned int RTree<R,I>::get_size_tree(){
  return get_size_tree(root);
}
template <typename R, typename I>
unsigned int RTree<R,I>::get_size_leaves(){
  return get_size_leaves(root);
}

template <typename R, typename I>
unsigned int RTree<R,I>::get_size_tree(const RTreeNode<R,I>* node)
{
  if (node == nullptr) return 0;
  unsigned int size = 1;
  if (node->entries.index() == 0) { // leaf node
    return size;
  } else {
    queue<const RTreeNode<R,I>*> q;
    q.push(node);

    while (q.size() > 0) {
      const RTreeNode<R,I>* next = q.front();
      q.pop();
      if (next->entries.index() == 0) { // next is a leaf
	size++;
	continue;
      }
      for (RTreeNode<R,I>* entry : std::get<1>(next->entries)) {
	size++;
	q.push(entry);
      }
    }
  }
  return size;
}

template <typename R, typename I>
unsigned int RTree<R,I>::get_size_leaves(const RTreeNode<R,I>* node)
{
  if (node == nullptr) return 0;
  unsigned int size = 1;
  if (node->entries.index() == 0) { // leaf node
    return size;
  } else {
    queue<const RTreeNode<R,I>*> q;
    q.push(node);

    while (q.size() > 0) {
      const RTreeNode<R,I>* next = q.front();
      q.pop();
      if (next->entries.index() == 0) { // next is a leaf
	size+=std::get<0>(next->entries).size();
	continue;
      }
      for (RTreeNode<R,I>* entry : std::get<1>(next->entries)) {
	q.push(entry);
      }
    }
  }
  return size;
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
  double max_d = -1e30;
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
  //std::cout << " got here afterall" << std::endl;
  //std::cout << seeds[0] << " " << seeds[1] << std::endl;
  std::set<unsigned int> g0 = {seeds[0]}, g1 = {seeds[1]};
  R g0mbr = *mbrs[seeds[0]], g1mbr = *mbrs[seeds[1]];

  int num_remaining = get_num_mbrs(node) - 2;

  // allocate all entries to either of the seeds preserving the minimum entries per node requirement
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
    node->mbr = g0mbr; // adjust node's mbr to reflect its current entries

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
void RTree<R,I>::insert(R mbr, I st_index)
{
  total_num_entries++;
  if (root == nullptr) {
    root = new RTreeNode<R,I>(mbr, nullptr, std::vector<LeafEntry<R,I>>({{mbr,st_index}}) );
    return;
  }

  RTreeNode<R,I>* target_leaf = choose_leaf(mbr);
  //std::cout << "	got leaf"<<std::endl;
  std::vector<LeafEntry<R,I>>& entries = std::get<0>(target_leaf->entries);
  entries.push_back({mbr, st_index});
  target_leaf->mbr = merge_f( target_leaf->mbr, mbr );

  RTreeNode<R,I>* twin_node = nullptr;
  if (entries.size() > max_entries) {
    //std::cout << "	splitting"<<std::endl;
    twin_node = split_node(target_leaf);
    //std::cout << "	splitted"<<std::endl;
  }
  //std::cout << "	adjusting"<<std::endl;
  adjust_tree(target_leaf);
  //std::cout << "	adjusted"<<std::endl;
    
}
/* *************************** R Tree Search methods ************************* */

template <typename R, typename I>
std::vector<I> RTree<R,I>::sim_search(const std::vector<double>& q, double epsilon)
{
  epsilon = epsilon * epsilon;
  if (uncompr_point_dist_sqr_f(q, root->mbr) > epsilon) {
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
	if (uncompr_point_dist_sqr_f(q, l.mbr) <= epsilon) {
	  results.push_back(l.st_index);
	}
      }
    } else {
      for (const RTreeNode<R,I>* const n : std::get<1>(next->entries)) {
	if (uncompr_point_dist_sqr_f(q, n->mbr) <= epsilon) {
	  nodes.push(n);
	}
      }
    }
  }
  return results;
}

template <typename R, typename I>
std::vector<std::array<const double*,2>> RTree<R,I>::knn_search(const std::vector<double>& q, unsigned int k, FPtrRetrievalMethod<I> retrieve_f, const std::vector<double>& s)
{
  typedef std::variant<const LeafEntry<R,I>*, const RTreeNode<R,I>*> QEntry;
  auto entry_error = [this,&q](const QEntry& a) {
    const R& ambr = a.index()==0 ? std::get<0>(a)->mbr : std::get<1>(a)->mbr;
    return uncompr_point_dist_sqr_f(q,ambr);
  };
  auto cmp_qentry = [this,&q,&entry_error](const QEntry& a, const QEntry& b) {
    return entry_error(a) > entry_error(b);
  };
  typedef std::tuple<QEntry, double> QDEntry;
  auto cmp = [](const QDEntry& a, const QDEntry& b) {
    return std::get<1>(a) > std::get<1>(b);
  };
  using std::vector, std::array;
  typedef array<const double*,2> Subseq;
  typedef std::tuple<Subseq,double> SubseqWithE;
  auto ptr_error = [this,&q](const Subseq& s) { return error_measures::se_between_ptrs(q.data(), q.data()+q.size()-1,s[0],s[1]); };
  auto leaf_cmp = [](const SubseqWithE& sa, const SubseqWithE& sb) {
    return std::get<1>(sa) > std::get<1>(sb);
  };

  std::priority_queue<QDEntry, std::vector<QDEntry>, decltype(cmp)> pri_q(cmp);
  std::priority_queue<SubseqWithE, std::vector<SubseqWithE>, decltype(leaf_cmp)> candidates(leaf_cmp);
  std::vector<array<const double*,2>> results;
  pri_q.push( { root, entry_error(root) } );

  while (pri_q.size() != 0) {
    auto [next, next_error] = pri_q.top();
    pri_q.pop();
    while (candidates.size()!=0 && std::get<1>(candidates.top()) <= next_error ) {
      results.push_back( std::get<0>(candidates.top()) );
      candidates.pop();
      if (results.size() >= k) {
	return results;
      }
    }

    if (next.index() == 0) { // next is an entry
      for (auto s : retrieve_f( std::get<0>(next)->st_index, s ) ) {
	candidates.push({ s, ptr_error(s) });
      }
    } else {
      const RTreeNode<R,I>* node = std::get<1>(next);
      if (node->entries.index() == 0) { // next is a leaf node
	for ( const LeafEntry<R,I>& l : std::get<0>(node->entries) ){
	  pri_q.push({ &l, entry_error(&l) } ); // potential optimisation if you can retrieve the original DRT and use dist_LB
	}
      } else {
	for ( const RTreeNode<R,I>* e : std::get<1>(node->entries) ){
	  pri_q.push( { e, entry_error(e) } );
	}
      }
    }
  }
  return results;
}

template <typename R, typename I>
std::vector<std::array<const double*,2>> RTree<R,I>::sim_search_exact(const std::vector<double>& q, double epsilon, FPtrRetrievalMethod<I> retrieve_f, const std::vector<double>& s)
{
  epsilon = epsilon * epsilon;
  if (uncompr_point_dist_sqr_f(q, root->mbr) > epsilon) {
    return {};
  }

  queue<const RTreeNode<R,I>*> nodes;
  nodes.push(root);
  std::vector<std::array<const double*, 2>> results;

  while (nodes.size() != 0) {
    const RTreeNode<R,I>* const next = nodes.front();
    nodes.pop();
    if (next->entries.index() == 0) {
      for (const LeafEntry<R,I>& l : std::get<0>(next->entries)) {
	if (uncompr_point_dist_sqr_f(q, l.mbr) <= epsilon) {
	  for (const auto& [s_ptr,e_ptr] : retrieve_f(l.st_index,s)) {
	    if (error_measures::se_between_ptrs(q.data(), q.data()+q.size()-1, s_ptr, e_ptr) <= epsilon)
	      results.push_back({s_ptr,e_ptr});
	  }
	}
      }
    } else {
      for (const RTreeNode<R,I>* const n : std::get<1>(next->entries)) {
	if (uncompr_point_dist_sqr_f(q, n->mbr) <= epsilon) {
	  nodes.push(n);
	}
      }
    }
  }
  return results;
}
template <typename R, typename I>
double RTree<R,I>::pruning_power(const std::vector<double>& q, FPtrRetrievalMethod<I> retrieve_f, const std::vector<double>& s)
{
  typedef std::variant<const LeafEntry<R,I>*, const RTreeNode<R,I>*> QEntry;
  typedef std::tuple<QEntry,double> QDEntry;
  auto entry_error = [this,&q](const QEntry& a) {
    const R& ambr = a.index()==0 ? std::get<0>(a)->mbr : std::get<1>(a)->mbr;
    return uncompr_point_dist_sqr_f(q,ambr);
  };
  auto cmp = [](const QDEntry& a, const QDEntry& b) {
    return std::get<1>(a) > std::get<1>(b);
  };
  using std::vector, std::array;
  auto ptr_error = [this,&q](const array<const double*,2>& s) { return error_measures::se_between_ptrs(q.data(), q.data()+q.size()-1,s[0],s[1]); };
  auto leaf_cmp = [this,&q,&ptr_error](const array<const double*,2>& sa, const array<const double*,2>& sb) {
    return ptr_error(sa) > ptr_error(sb);
  };

  std::priority_queue<QDEntry, std::vector<QDEntry>, decltype(cmp)> pri_q(cmp);
  pri_q.push( { root, entry_error(root) } );

  array<const double*,2> min_index;
  double min_actual_error = 100000000000000000;
  double objects_scanned = 0.0;

  while (pri_q.size() != 0) {
    auto [next,next_error] = pri_q.top();
    pri_q.pop();
    if (min_actual_error < next_error) {
      return objects_scanned / (double) total_num_entries;
    }

    if (next.index() == 0) { // next is an entry
      objects_scanned+=1.0;
      for (auto s : retrieve_f( std::get<0>(next)->st_index, s ) ) {
	if (ptr_error(s) < min_actual_error) {
	  min_actual_error = ptr_error(s);
	  min_index = s;
	}
      }
    } else {
      const RTreeNode<R,I>* node = std::get<1>(next);
      if (node->entries.index() == 0) { // next is a leaf node
	for ( const LeafEntry<R,I>& l : std::get<0>(node->entries) ){
	  pri_q.push({ &l, entry_error(&l) }); // potential optimisation if you can retrieve the original DRT and use dist_LB
	}
      } else {
	for ( const RTreeNode<R,I>* e : std::get<1>(node->entries) ){
	  pri_q.push({e,entry_error(e)});
	}
      }
    }
  }
  return 1.0;
}

#endif
