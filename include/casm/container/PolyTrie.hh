
#ifndef POLYTRIE_HH
#define POLYTRIE_HH

#include <iostream>
#include <cassert>
#include <new>
#include <stdlib.h>

#include "casm/container/Array.hh"
#include "casm/misc/CASM_math.hh"

namespace CASM {

  template<typename T>
  class PolyTrie;
  template<typename T>
  class PTNode;
  template<typename T>
  class PTLeaf;


  template<typename T>
  class PTNode : protected Array<PTNode<T>* > {
  protected:
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
    friend class PolyTrie<T>;

    PTLeaf<T> *next_leaf;
    PTLeaf<T> **prev_leaf_addr;
    Array<Index> m_key;

  public:
    PTLeaf(PTNode<T> *_up, const Array<Index> &_key, const T &_val) :
      PTNode<T>(_up, _val), next_leaf(nullptr), prev_leaf_addr(nullptr), m_key(_key) {};

    ~PTLeaf();

    using PTNode<T>::remove;

    ///virtual from PTNode<T>::remove
    ///removes from list
    //void remove();

    const Array<Index> &key() const {
      return m_key;
    };

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
    void insert_after(PTLeaf<T> **prev);
    void insert_before(PTLeaf<T> *next);

  };

  //==========================================================

  template<typename T>
  class PolyTrie : public PTNode<T> {
    Index m_depth;
    PTLeaf<T> *m_leaf_list;
    using PTNode<T>::size;
    using PTNode<T>::at;
  public:

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

    /// set() allows assignment that prunes the trie if _val is approximately 0;
    void set(const Array<Index> &ind, const T &_val);

    void list_leaf(PTLeaf<T> *new_leaf);

    /// removes zero entries, if there are any, and returns true if entries were removed.
    bool prune_zeros(double tol = TOL);

    PTLeaf<T> *begin() {
      return m_leaf_list;
    };
    PTLeaf<T> const *begin() const {
      return m_leaf_list;
    };

    void print_sparse(std::ostream &out) const;

    //Arithmetic operations
    PolyTrie &operator*=(const T &scale);
    PolyTrie &operator+=(const PolyTrie<T> &RHS);
    PolyTrie &operator-=(const PolyTrie<T> &RHS);

    PolyTrie operator+(const PolyTrie<T> &RHS);
    PolyTrie operator-(const PolyTrie<T> &RHS);
  };

  //==========================================================

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
      if(!home_trie.begin())
        home_trie.list_leaf(new PTLeaf<T>(this, ind, 0));
      return home_trie.begin();
    };

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
  void PTLeaf<T>::insert_after(PTLeaf<T> **ptr_in_prev) {
    // 'ptr_in_prev' points to 'next_leaf' in the PTLeaf that will become 'previous'
    // This is the definition of prev_leaf_addr
    prev_leaf_addr = ptr_in_prev;

    // (*ptr_in_prev) currently points to leaf that will become 'next'
    next_leaf = *ptr_in_prev;

    // next_leaf->prev_leaf_addr should point to the next_leaf pointer in this leaf
    if(next_leaf)
      next_leaf->prev_leaf_addr = &next_leaf;

    // (*ptr_in_prev) is next_leaf in 'previous', and it should be set to this
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
    PTLeaf<T> const *current(orig.begin());
    //PolyTrie<T>::set() does all the work -- result is a pruned PolyTrie that is a copy of 'orig'
    while(current) {
      set(current->key(), current->val());
      current = current->next();
    }

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~PTNode() should take care of PolyTrie destruction, except for the edge case m_depth=0
  template<typename T>
  PolyTrie<T>::~PolyTrie() {
    if(depth() == 0 && begin()) {
      delete begin();
    }

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  T PolyTrie<T>::get(const Array<Index> &ind) const {
    assert(ind.size() == depth() && "In PolyTrie<T>::get(), ind.size() must match PolyTrie<T>::depth().");
    if(ind.size() == 0) {
      if(begin())
        return begin()->val();
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
    assert(ind.size() == depth() && "In PolyTrie<T>::get(), ind.size() must match PolyTrie<T>::depth().");

    if(almost_zero(_val)) {
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
    assert(ind.size() == depth() && "In PolyTrie<T>::get(), ind.size() must match PolyTrie<T>::depth().");

    PTNode<T> *tnode(this);
    for(Index i = 0; i < ind.size() - 1; i++) {
      tnode = tnode->valid_node_at(ind[i]);
    }

    tnode = tnode->valid_leaf_at(*this, ind);

    return tnode->m_val;

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  void PolyTrie<T>::list_leaf(PTLeaf<T> *new_leaf) {
    new_leaf->insert_after(&m_leaf_list);
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  void PolyTrie<T>::clear() {
    PTLeaf<T> *current(begin());
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

    std::swap(m_leaf_list, other.m_leaf_list);
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
    //std::cout << "Entering prune_zeros!!!!!!!!!!!!!!!!\n";
    bool is_edited(false);
    PTLeaf<T> *current(begin());
    while(current) {
      if(almost_zero(current->val(), tol)) {
        current = current->remove_and_next();
        is_edited = true;
      }
      else {
        current = current->next();
      }
    }
    return is_edited;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  void PolyTrie<T>::print_sparse(std::ostream &out) const {
    PTLeaf<T> const *current(begin());

    while(current) {
      if(!almost_zero(current->val())) {
        out << current-> key() << ":  " << current->val() << '\n';
      }
      current = current->next();

    }


  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  PolyTrie<T> &PolyTrie<T>::operator*=(const T &scale) {
    if(almost_zero(scale)) {
      clear();
      return *this;
    }

    PTLeaf<T> *current(begin());
    while(current) {
      (current->m_val) *= scale;
      current = current->next();
    }
    return *this;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  PolyTrie<T> &PolyTrie<T>::operator+=(const PolyTrie<T> &RHS) {
    PTLeaf<T> const *current(RHS.begin());
    while(current) {
      if(almost_zero(at(current->key()) += current->val())) {
        remove(current->key());
      }
      current = current->next();
    }
    return *this;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  PolyTrie<T> &PolyTrie<T>::operator-=(const PolyTrie<T> &RHS) {
    PTLeaf<T> const *current(RHS.begin());
    while(current) {
      if(almost_zero(at(current->key()) -= current->val())) {
        remove(current->key());
      }
      current = current->next();
    }
    return *this;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  PolyTrie<T> PolyTrie<T>::operator+(const PolyTrie<T> &RHS) {
    PolyTrie<T> result(*this);
    return result += RHS;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  PolyTrie<T> PolyTrie<T>::operator-(const PolyTrie<T> &RHS) {
    PolyTrie<T> result(*this);
    return result -= RHS;
  }
}

#endif /* POLYTRIE_H_ */
