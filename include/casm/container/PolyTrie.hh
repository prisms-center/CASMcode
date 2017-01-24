
#ifndef POLYTRIE_HH
#define POLYTRIE_HH

#include <iostream>
#include <cassert>
#include <new>
#include <stdlib.h>

#include "casm/CASM_global_definitions.hh"
#include "casm/container/Array.hh"
#include "casm/misc/CASM_TMP.hh"

namespace CASM {

  /**
   *  \defgroup PolyTrie
   *
   *  \ingroup Container
   *
   *  \brief A trie structure for representing polynomials
   *
   *  @{
   */

  template<typename T>
  class PolyTrie;
  template<typename T>
  class PTNode;
  template<typename T>
  class PTLeaf;

  namespace ComparePTLeaf {
    class ByMonomialOrder;
  }
  template<typename PTType, bool IsConst = false>
  class PTIterator;

  template<typename T>
  class PTNode : protected Array<PTNode<T>* > {
  protected:
    static double PT_TOL() {
      return 1e-6;
    }
    friend class PolyTrie<T>;
    PTNode<T> *up_node;
    T m_val;

    using Array<PTNode<T>* > :: size;
    using Array<PTNode<T>* > :: at;
    using Array<PTNode<T>* > :: push_back;

    PTNode<T> *valid_node_at(Index i);
    PTNode<T> *valid_leaf_at(PolyTrie<T> &home_trie, const Array<Index> &ind);

  public:

    PTNode(PTNode<T> *_up, const T &_val = 0) :
      up_node(_up),  m_val(_val) {};

    const T &val()const {
      return m_val;
    };
    virtual ~PTNode();

    virtual void remove();
  };

  //==========================================================
  template<typename T>
  class PTLeaf : public PTNode<T> {
  private:
    PTLeaf<T> *next_leaf;
    PTLeaf<T> **prev_leaf_addr;
    Array<Index> m_key;

  protected:
    friend class PolyTrie<T>;

  public:
    PTLeaf(PTNode<T> *_up, const Array<Index> &_key, const T &_val) :
      PTNode<T>(_up, _val), next_leaf(nullptr), prev_leaf_addr(nullptr), m_key(_key) {};

    ~PTLeaf();

    typedef ComparePTLeaf::ByMonomialOrder CompareByMonomialOrder;

    using PTNode<T>::remove;

    ///virtual from PTNode<T>::remove
    ///removes from list
    //void remove();

    const Array<Index> &key() const {
      return m_key;
    };

    T &val_ref() {
      return PTNode<T>::m_val;
    }
    const T &val_ref() const {
      return PTNode<T>::m_val;
    }

    //Const and non-const list traversal
    PTLeaf<T> *next() {
      return next_leaf;
    };
    PTLeaf<T> const *next()const {
      return next_leaf;
    };

    PTLeaf<T> *remove_and_next();

    /// Do a swap with 'other' and 'this' so that leaf previous to 'other' now points to 'this'
    /// and leaf previous to 'this' now points to 'other'
    void swap_after_prev(PTLeaf<T> *other);

    ///List insertion methods
    //insert at is passed a reference to the pointer at which insertion will occur,
    // e.g.: ptr_to_new_leaf->insert_at(ptr_to_list_head); <-- ptr_to_list_head can be NULL
    void insert_at(PTLeaf<T> *&insertion_ptr);

    // prev can't be NULL
    void insert_after(PTLeaf<T> *prev);
    // next can't be NULL
    void insert_before(PTLeaf<T> *next);
    void make_tail() {
      next_leaf = NULL;
    }

  };

  namespace ComparePTLeaf {
    class ByMonomialOrder {
    public:
      template <typename T>
      static bool compare(const PTLeaf<T> &A, const PTLeaf<T> &B) {
        return (A.key().sum() >= B.key().sum()) && (A.key() >= B.key());
      }
      template <typename T>
      bool operator()(const PTLeaf<T> &A, const PTLeaf<T> &B) const {
        return compare(A, B);
      }
    };

    class CustomOrder {
    public:
      template <typename T>
      static bool compare(const PTLeaf<T> &A, const PTLeaf<T> &B) {
        if(A.key().size() < 9 || B.key().size() < 9)
          return ByMonomialOrder::compare(A, B);

        if(A.key().sum() > B.key().sum())
          return true;
        if(A.key().sum() < B.key().sum())
          return false;

        Array<Index>
        sub1A(A.key().sub_array(0, A.key().size() - 7)),
              sub2A(A.key().sub_array(A.key().size() - 6, A.key().size() - 4)),
              sub3A(A.key().sub_array(A.key().size() - 3, A.key().size() - 1)),
              sub1B(B.key().sub_array(0, B.key().size() - 7)),
              sub2B(B.key().sub_array(B.key().size() - 6, B.key().size() - 4)),
              sub3B(B.key().sub_array(B.key().size() - 3, B.key().size() - 1));
        if(sub1A.sum() > sub1B.sum())
          return true;
        if(sub1A.sum() < sub1B.sum())
          return false;
        if(sub1A.max() > sub1B.max())
          return true;
        if(sub1A.max() < sub1B.max())
          return false;
        if(sub2A.sum() > sub2B.sum())
          return true;
        if(sub2A.sum() < sub2B.sum())
          return false;
        if(sub3A.sum() > sub3B.sum())
          return true;
        if(sub3A.sum() < sub3B.sum())
          return false;
        if(almost_equal(std::abs(A.val()), std::abs(B.val())))
          return ByMonomialOrder::compare(A, B) && A.val() + TOL > B.val();

        return std::abs(A.val()) > std::abs(B.val());
      }
      template <typename T>
      bool operator()(const PTLeaf<T> &A, const PTLeaf<T> &B) const {
        return compare(A, B);
      }
    };
  }

  //==========================================================

  template<typename T>
  class PolyTrie : public PTNode<T> {
    Index m_depth;
    PTLeaf<T> *m_leaf_list;
    using PTNode<T>::size;
    using PTNode<T>::at;

    Index num_nonzero() const;

    using PTNode<T>::valid_node_at;
    using PTNode<T>::valid_leaf_at;
  public:
    typedef T value_type;

    typedef PTLeaf<T> leaf_type;

    typedef PTIterator<PolyTrie<T> > iterator;
    typedef PTIterator<PolyTrie<T>, true> const_iterator;

    PolyTrie(Index _depth) : PTNode<T>(nullptr), m_depth(_depth), m_leaf_list(nullptr) {};
    PolyTrie(const PolyTrie<T> &orig);
    // ~PTNode() should take care of PolyTrie destruction, except for the edge case m_depth=0
    ~PolyTrie();

    ///Virtual from PTNode<T>, you cannot remove a PolyTrie
    void remove() { };

    ///Remove entry at 'ind'
    void remove(const Array<Index> &ind);

    /// Efficient swap of this PolyTrie with another
    void swap(PolyTrie<T> &other);

    ///Clears all leaves and resets depth
    void redefine(Index new_depth);

    void clear();

    Index depth() const {
      return m_depth;
    };

    ///get() provides constant access does not change trie structure
    T get(const Array<Index> &ind)const;

    /// at() provides non-const access and changes trie if leaf does not exist at ind
    T &at(const Array<Index> &ind);

    T &operator()(const Array<Index> &ind);

    T operator()(const Array<Index> &ind) const;

    /// set() allows assignment that prunes the trie if _val is approximately 0;
    void set(const Array<Index> &ind, const T &_val);

    void list_leaf(PTLeaf<T> *new_leaf);

    /// removes zero entries, if there are any, and returns true if entries were removed.
    bool prune_zeros(double tol = PTNode<T>::PT_TOL());

    PTLeaf<T> *begin_ptr() {
      return m_leaf_list;
    }
    PTLeaf<T> const *begin_ptr() const {
      return m_leaf_list;
    }

    iterator begin() {
      return iterator(m_leaf_list);
    }
    const_iterator begin() const {
      return const_iterator(m_leaf_list);
    }

    iterator end() {
      return iterator();
    }
    const_iterator end() const {
      return const_iterator();
    };

    void print_sparse(std::ostream &out) const;

    //Arithmetic operations
    PolyTrie &operator*=(const T &scale);
    PolyTrie &operator+=(const PolyTrie<T> &RHS);
    PolyTrie &operator-=(const PolyTrie<T> &RHS);

    PolyTrie operator+(const PolyTrie<T> &RHS) const;
    PolyTrie operator-(const PolyTrie<T> &RHS) const;

    bool compare(const PolyTrie<T> &RHS, double tol = 2 * PTNode<T>::PT_TOL()) const;
    bool operator==(const PolyTrie<T> &RHS) const;
    bool almost_zero(double tol = 2 * PTNode<T>::PT_TOL()) const;

    template<typename CompareType>
    void sort_leaves(const CompareType &compare);
  };

  //==========================================================
  template<typename PTType, bool IsConst>
  class PTIterator {
  public:
    typedef typename CASM_TMP::ConstSwitch<IsConst, typename PTType::leaf_type> leaf_type;
    typedef typename CASM_TMP::ConstSwitch<IsConst, typename PTType::value_type> *pointer;
    typedef typename CASM_TMP::ConstSwitch<IsConst, typename PTType::value_type> &reference;

  private:
    leaf_type *m_curr_leaf_ptr;

  public:
    PTIterator(leaf_type *_leaf_ptr = 0) : m_curr_leaf_ptr(_leaf_ptr) {};
    PTIterator(const PTIterator<PTType, false> &_it);

    template<bool IsConst2>
    bool operator==(const PTIterator<PTType, IsConst2> &_it) {
      return m_curr_leaf_ptr == _it.m_curr_leaf_ptr;
    }

    template<bool IsConst2>
    bool operator!=(const PTIterator<PTType, IsConst2> &_it) {
      return !((*this) == _it);
    }

    PTIterator &operator++() {
      if(m_curr_leaf_ptr)
        m_curr_leaf_ptr = m_curr_leaf_ptr->next();
      return *this;
    }
    PTIterator operator++(int) {
      PTIterator t_it(*this);
      ++(*this);
      return t_it;
    }

    PTIterator &remove_and_next() {
      m_curr_leaf_ptr = m_curr_leaf_ptr->remove_and_next();
      return *this;
    }

    reference operator*() {
      return m_curr_leaf_ptr->val_ref();
    }
    pointer operator->() {
      return &operator*;
    }
    const Array<Index> &key() {
      return m_curr_leaf_ptr->key();
    }

    leaf_type *leaf_ptr() const {
      return m_curr_leaf_ptr;
    };

  };

  //==========================================================

  template<typename PTType, bool IsConst>
  PTIterator<PTType, IsConst>::PTIterator(const PTIterator<PTType, false> &_it) :
    m_curr_leaf_ptr(_it.leaf_ptr()) {}

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  PTNode<T> *PTNode<T>::valid_node_at(Index i) {
    while(size() <= i) {
      push_back(nullptr);
    }
    return at(i) ? at(i) : at(i) = new PTNode<T>(this);
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  PTNode<T> *PTNode<T>::valid_leaf_at(PolyTrie<T> &home_trie, const Array<Index> &ind) {
    if(ind.size() == 0 && home_trie.depth() == 0 && this == static_cast<PTNode<T> *>(&home_trie)) {
      if(!home_trie.begin_ptr())
        home_trie.list_leaf(new PTLeaf<T>(this, ind, 0));
      return home_trie.begin_ptr();
    }

    while(size() <= ind.back()) {
      push_back(nullptr);
    }

    // Check if there's a valid leaf at 'ind' element and create a new one if not
    if(!at(ind.back())) {
      at(ind.back()) = new PTLeaf<T>(this, ind, 0);
      home_trie.list_leaf(static_cast<PTLeaf<T>*>(at(ind.back())));
    }

    return at(ind.back());
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  void PTNode<T>::remove() {

    if(!up_node) {
      delete this;
      return;
    }

    bool remove_up(true), not_removed(true);

    for(Index i = 0; i < up_node->size() && (remove_up || not_removed); i++) {
      if((up_node->at(i)) == this) {
        up_node->at(i) = nullptr;
        not_removed = false;
      }
      remove_up = remove_up && !(up_node->at(i));
    }
    if(remove_up) {
      //clearing up_node saves some time when it is deleted
      up_node->clear();
      up_node->remove();
    }
    //std::cout << "Finishing remove by deleting myself:\n";
    delete this;
    return;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  PTNode<T>::~PTNode<T>() {
    for(Index i = 0; i < size(); i++) {
      if(at(i))
        delete at(i);
    }

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  void PTLeaf<T>::swap_after_prev(PTLeaf<T> *other) {
    swap(prev_leaf_addr, other->prev_leaf_addr);
    *prev_leaf_addr = this;
    *(other->prev_leaf_addr) = other;

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  template<typename T>
  void PTLeaf<T>::insert_after(PTLeaf<T> *prev) {
    insert_at(prev->next_leaf);
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  void PTLeaf<T>::insert_at(PTLeaf<T> *&insertion_ptr) {
    // 'insertion_ptr' refers to 'next_leaf' in the PTLeaf that will become 'previous' (or pointer that is head of list)
    // This is the definition of prev_leaf_addr
    prev_leaf_addr = &insertion_ptr;

    // insertion_ptr currently points to leaf that will become 'next'
    next_leaf = insertion_ptr;

    // next_leaf->prev_leaf_addr should point to the next_leaf pointer in this leaf
    if(next_leaf)
      next_leaf->prev_leaf_addr = &next_leaf;

    // insertion_ptr is next_leaf in 'previous', and it should be set to this
    // the following line is equivalent to insertion_ptr=this
    (*prev_leaf_addr) = this;


  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  void PTLeaf<T>::insert_before(PTLeaf<T> *next) {
    assert(next && next->prev_leaf_addr && "PTLeaf insertion failed because 'next' leaf address is not in List.");

    // 'next' will become 'next_leaf'
    next_leaf = next;

    // prev_leaf_addr should be same as next->prev_leaf_addr
    prev_leaf_addr = next->prev_leaf_addr;

    // *prev_leaf_addr is next_leaf in 'previous' leaf; it should point to this
    (*prev_leaf_addr) = this;

    // next->prev_leaf_addr should point to next_leaf ptr in this
    next->prev_leaf_addr = &next_leaf;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  template<typename T>
  PTLeaf<T>::~PTLeaf<T>() {
    //std::cout << "I'm PTLeaf " << this << " and I'm being destroyed!!! next_leaf is " << next_leaf << " and prev_leaf_addr is " << prev_leaf_addr << " and " << *prev_leaf_addr << "\n";

    if(prev_leaf_addr)
      *prev_leaf_addr = next_leaf;

    if(next_leaf) {
      next_leaf->prev_leaf_addr = prev_leaf_addr;
      next_leaf = nullptr;
    }
    prev_leaf_addr = nullptr;

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  template<typename T>
  PTLeaf<T> *PTLeaf<T>::remove_and_next() {
    PTLeaf<T> *save_next(next_leaf);

    //remove deletes 'this' so can only return afterwards
    remove();
    return save_next;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  PolyTrie<T>::PolyTrie(const PolyTrie<T> &orig) : PTNode<T>(nullptr), m_depth(orig.m_depth), m_leaf_list(nullptr) {
    PolyTrie<T>::const_iterator it(orig.begin()), it_end(orig.end());
    //PolyTrie<T>::set() does all the work -- result is a pruned PolyTrie that is a copy of 'orig'
    for(; it != it_end; ++it)
      set(it.key(), *it);
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~PTNode() should take care of PolyTrie destruction, except for the edge case m_depth=0
  template<typename T>
  PolyTrie<T>::~PolyTrie() {
    if(depth() == 0 && begin_ptr())
      delete begin_ptr();
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  Index PolyTrie<T>::num_nonzero() const {
    Index count(0);
    PolyTrie<T>::const_iterator it(begin()), it_end(end());
    for(; it != it_end; ++it)
      count++;
    return count;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  T PolyTrie<T>::get(const Array<Index> &ind) const {
    assert(ind.size() == depth() && "In PolyTrie<T>::get(), ind.size() must match PolyTrie<T>::depth().");
    if(ind.size() == 0) {
      if(begin_ptr())
        return begin_ptr()->val();
      return 0;
    }

    PTNode<T> const *tnode(this);
    for(Index i = 0; i < ind.size(); i++) {
      if(ind[i] >= tnode->size() || !(tnode->at(ind[i]))) {
        return 0;
      }
      tnode = tnode->at(ind[i]);
    }
    return tnode->val();
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  void PolyTrie<T>::set(const Array<Index> &ind, const T &_val) {
    assert(ind.size() == depth() && "In PolyTrie<T>::set(), ind.size() must match PolyTrie<T>::depth().");

    if(CASM::almost_zero(_val, PTNode<T>::PT_TOL())) {
      remove(ind);
    }

    PTNode<T> *tnode(this);
    for(Index i = 0; i < ind.size() - 1; i++) {
      tnode = tnode->valid_node_at(ind[i]);
    }

    tnode = tnode->valid_leaf_at(*this, ind);
    tnode->m_val = _val;
    return;

  }


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  T &PolyTrie<T>::at(const Array<Index> &ind) {
    assert(ind.size() == depth() && "In PolyTrie<T>::at(), ind.size() must match PolyTrie<T>::depth().");

    PTNode<T> *tnode(this);
    for(Index i = 0; i < ind.size() - 1; i++) {
      tnode = tnode->valid_node_at(ind[i]);
    }

    tnode = tnode->valid_leaf_at(*this, ind);

    return tnode->m_val;

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  template<typename T>
  T &PolyTrie<T>::operator()(const Array<Index> &ind) {
    assert(ind.size() == depth() && "In PolyTrie<T>::operator(), ind.size() must match PolyTrie<T>::depth().");

    PTNode<T> *tnode(this);
    for(Index i = 0; i < ind.size() - 1; i++) {
      tnode = tnode->valid_node_at(ind[i]);
    }

    tnode = tnode->valid_leaf_at(*this, ind);

    return tnode->m_val;

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  template<typename T>
  T PolyTrie<T>::operator()(const Array<Index> &ind) const {
    assert(ind.size() == depth() && "In PolyTrie<T>::get(), ind.size() must match PolyTrie<T>::depth().");
    if(ind.size() == 0) {
      if(begin_ptr())
        return begin_ptr()->val();
      return 0;
    }

    PTNode<T> const *tnode(this);
    for(Index i = 0; i < ind.size(); i++) {
      if(ind[i] >= tnode->size() || !(tnode->at(ind[i]))) {
        return 0;
      }
      tnode = tnode->at(ind[i]);
    }
    return tnode->val();


  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  void PolyTrie<T>::list_leaf(PTLeaf<T> *new_leaf) {
    new_leaf->insert_at(m_leaf_list);
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  void PolyTrie<T>::clear() {
    PTLeaf<T> *current(begin_ptr());
    while(current) {
      current = current->remove_and_next();
    }
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  void PolyTrie<T>::redefine(Index new_depth) {
    clear();
    m_depth = new_depth;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  void PolyTrie<T>::swap(PolyTrie<T> &other) {
    if(this == &other)
      return;

    for(Index i = 0; i < size(); i++) {
      if(at(i)) {
        at(i)->up_node = &other;
      }
    }
    Array<PTNode<T>*>::swap(other);

    for(Index i = 0; i < size(); i++) {
      if(at(i)) {
        at(i)->up_node = this;
      }
    }

    CASM::swap(m_leaf_list, other.m_leaf_list);
    if(m_leaf_list)
      m_leaf_list->prev_leaf_addr = &m_leaf_list;
    if(other.m_leaf_list)
      other.m_leaf_list->prev_leaf_addr = &(other.m_leaf_list);

  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  void PolyTrie<T>::remove(const Array<Index> &ind) {
    assert(ind.size() == depth() && "In PolyTrie<T>::remove(), ind.size() must match PolyTrie<T>::depth().");
    PTNode<T> *tnode(this);
    for(Index i = 0; i < ind.size(); i++) {
      if(ind[i] >= tnode->size() || !(tnode->at(ind[i]))) {
        return;
      }
      tnode = tnode->at(ind[i]);
    }
    tnode->remove();
    return;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  bool PolyTrie<T>::prune_zeros(double tol) {
    bool is_edited(false);
    PolyTrie<T>::iterator it(begin()), it_end(end());
    for(; it != it_end; ++it) {
      while(it != it_end && CASM::almost_zero(*it, tol)) {
        it.remove_and_next();
        is_edited = true;
      }
    }
    return is_edited;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  void PolyTrie<T>::print_sparse(std::ostream &out) const {
    PolyTrie<T>::const_iterator it(begin()), it_end(end());
    for(; it != it_end; ++it) {
      if(!CASM::almost_zero(*it, PTNode<T>::PT_TOL()))
        out << it.key() << ":  " << *it << '\n';
    }
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  PolyTrie<T> &PolyTrie<T>::operator*=(const T &scale) {
    if(CASM::almost_zero(scale, PTNode<T>::PT_TOL())) {
      clear();
      return *this;
    }

    PolyTrie<T>::iterator it(begin()), it_end(end());
    for(; it != it_end; ++it)
      *it *= scale;

    return *this;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  PolyTrie<T> &PolyTrie<T>::operator+=(const PolyTrie<T> &RHS) {
    PolyTrie<T>::const_iterator it(RHS.begin()), it_end(RHS.end());
    for(; it != it_end; ++it) {
      if(CASM::almost_zero(at(it.key()) += *it, PTNode<T>::PT_TOL()))
        remove(it.key());
    }
    return *this;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  PolyTrie<T> &PolyTrie<T>::operator-=(const PolyTrie<T> &RHS) {
    PolyTrie<T>::const_iterator it(RHS.begin()), it_end(RHS.end());
    for(; it != it_end; ++it) {
      if(CASM::almost_zero(at(it.key()) -= *it, PTNode<T>::PT_TOL()))
        remove(it.key());
    }
    return *this;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  PolyTrie<T> PolyTrie<T>::operator+(const PolyTrie<T> &RHS) const {
    PolyTrie<T> result(*this);
    return result += RHS;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  PolyTrie<T> PolyTrie<T>::operator-(const PolyTrie<T> &RHS) const {
    PolyTrie<T> result(*this);
    return result -= RHS;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  bool PolyTrie<T>::compare(const PolyTrie<T> &RHS, double tol) const {
    return ((*this) - RHS).almost_zero(tol);
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  bool PolyTrie<T>::operator==(const PolyTrie<T> &RHS) const {
    return compare(RHS);
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  bool PolyTrie<T>::almost_zero(double tol) const {
    PolyTrie<T>::const_iterator it(begin()), it_end(end());

    for(; it != it_end; ++it) {
      if(!CASM::almost_zero(*it, tol))
        return false;
    }
    return true;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T> template<typename CompareType>
  void PolyTrie<T>::sort_leaves(const CompareType &compare) {
    if(!m_leaf_list)
      return;

    //Working with two sublists, 'a_list', and 'b_list' -- define pointers to each
    PTLeaf<T> *a_ptr, *b_ptr;

    // merge_elem is the element to be merged and tail is the merge position
    PTLeaf<T> *merge_elem, *tail;

    // size of a_list and size of b_list
    int a_size, b_size;

    // merge_size is the target size of a_list and b_list at the beggining of each iteration.
    //nmerges is the number of merges performed during the interation
    int merge_size(1), nmerges(2);
    // count number of merges we do in this pass
    // If we have done only one merge, we're finished.
    while(nmerges > 1) {
      nmerges = 0;

      a_ptr = m_leaf_list;
      tail = NULL;

      while(a_ptr) {
        // there exists a merge to be done
        nmerges++;
        // step `merge_size' places along from p
        b_ptr = a_ptr;
        a_size = 0;
        for(int i = 0; i < merge_size; i++) {
          a_size++;
          if(!(b_ptr = b_ptr->next())) break;
        }

        // if b_ptr hasn't fallen off the end of the list, we have two m_leaf_lists to merge
        b_size = merge_size;

        // now we merge a_list with b_list, keeping the result ordered
        while(a_size > 0 || (b_size > 0 && b_ptr)) {
          // decide whether next merge_elem comes from a_list or b_list
          if(a_size == 0) {
            // a_list is empty; merge_elem must come from b_list.
            merge_elem = b_ptr;
            b_ptr = b_ptr->next();
            b_size--;
          }
          else if(b_size == 0 || !b_ptr) {
            // b_list is empty; merge_elem must come from a_list.
            merge_elem = a_ptr;
            a_ptr = a_ptr->next();
            a_size--;
          }
          else if(compare(*a_ptr, *b_ptr)) {
            // First element of b_list is "truer" (or same); merge_elem must come from a_list.
            merge_elem = a_ptr;
            a_ptr = a_ptr->next();
            a_size--;
          }
          else {
            // First element of b_list is "truer"; merge_elem must come from b_list.
            merge_elem = b_ptr;
            b_ptr = b_ptr->next();
            b_size--;
          }

          // add the next PTLeaf<T> to the merged list
          if(tail) {
            merge_elem->insert_after(tail);
          }
          else {
            merge_elem->insert_at(m_leaf_list);
          }
          tail = merge_elem;
        }

        // we've reached the end of a_list and b_list, so point a_ptr at the the same position as b_ptr,
        // which will either be NULL or the beginning of the next a_list
        a_ptr = b_ptr;
      }
      // set the new element as the tail (set internal pointer to NULL)
      tail->make_tail();

      //double the merge size for next iteration
      merge_size *= 2;
    }


  }

  /** @} */

}

#endif /* POLYTRIE_H_ */
