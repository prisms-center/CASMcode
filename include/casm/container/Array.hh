
#ifndef ARRAY_HH
#define ARRAY_HH

#include <iostream>
#include <cassert>
#include <new>
#include <stdlib.h>

#include "casm/CASM_global_definitions.hh"
#include "casm/misc/CASM_TMP.hh"
#include "casm/casm_io/jsonParser.hh"

namespace CASM {
  template<class T>
  class Array;

  template<class T>
  class ReturnArray;

  /**
   *  \defgroup Container
   *
   *  \brief Useful containers
   *
   */

  /**
   *  \defgroup Array
   *
   *  \ingroup Container
   *  \brief Basic std::vector like container (deprecated)
   *
   *  @{
   */

  /// \brief Basic std::vector like container (deprecated)
  ///
  /// - Legacy container that incorporated various useful algorithms, now being
  ///   slowly replaced
  /// - Prefer to use std::vector except when necessary for compatibility
  ///
  template<class T>
  class Array {
  private:
    static Index ARRAY_MIN_EXTRA_SPACE() {
      return 2;
    }
    static double ARRAY_EXTENSION_FACTOR() {
      return 1.1;
    }

    //static const int min_extra_space = 2;
    //static const double percent_extra_space = 0.1;


    Index N;
    Index NMax;
    T *Vals;

  public:
    //typedefs for nested Arrays -- Declare using 'Array<T>::XN', where 'N' is the number of nested indices
    typedef Array<T> X1;
    typedef Array<X1> X2;
    typedef Array<X2> X3;
    typedef Array<X3> X4;
    typedef Array<X4> X5;
    typedef Array<X5> X6;
    typedef Array<X6> X7;
    typedef Array<X7> X8;
    typedef Array<X8> X9;

    typedef T value_type;
    typedef Index size_type;
    typedef T *iterator;
    typedef const T *const_iterator;

    // CONSTRUCTORS/DESTRUCTORS
    Array() : N(0), NMax(0), Vals(nullptr) { }

    //*******************************************************************************************

    Array(Index init_N) : N(0), NMax(0), Vals(nullptr) {
      resize(init_N);
    }

    //*******************************************************************************************

    Array(Index init_N, const T &init_val) : N(0), NMax(0), Vals(nullptr) {
      resize(init_N, init_val);
    }

    //*******************************************************************************************

    Array(const Array &RHS) : N(0), NMax(0), Vals(nullptr) {
      reserve(RHS.size());
      for(Index i = 0; i < RHS.size(); i++)
        push_back(RHS[i]);
    }

    //*******************************************************************************************

    template<typename Iterator>
    Array(Iterator begin,
          Iterator end,
          typename CASM_TMP::enable_if_iterator<Iterator>::type * = nullptr) :
      N(0), NMax(0), Vals(nullptr) {

      reserve(std::distance(begin, end));
      auto it = begin;
      for(; it != end; ++it)
        push_back(*it);
    }

    //*******************************************************************************************
    Array(std::initializer_list<T> in) :
      N(0), NMax(0), Vals(nullptr) {

      reserve(in.size());
      auto it = in.begin();
      for(; it != in.end(); ++it)
        push_back(*it);
    }

    //*******************************************************************************************

    Array(ReturnArray<T> &RHS);


    ~Array() {
      clear();
      if(Vals)
        operator delete(Vals);
    }

    /// Returns an array with the sequence (initial, ++initial, ..., final), inclusive
    /// requires that operator<() and operator++() are defined on type T
    static ReturnArray<T> sequence(const T &initial, const T &final);

    /// Returns an array with the sequence (initial, initial+increment, ..., final?),
    /// inclusive if final is in the sequence
    /// requires that operator<() and operator+=() are defined on type T
    static ReturnArray<T> sequence(const T &initial, const T &increment, const T &final);

    Index size() const {
      return N;
    }

    // ASSIGN/REASSIGN
    Array &operator=(const Array &RHS);
    Array &operator=(ReturnArray<T> &RHS);
    void swap(Array<T> &RHS);


    // ACCESSORS

    T &at(Index ind) {
      assert(ind >= 0 && ind < N);
      return Vals[ind];
    }

    const T &at(Index ind) const {
      assert(ind >= 0 && ind < N);
      return Vals[ind];
    }

    const T &operator[](Index ind) const {
      assert(ind >= 0 && ind < N);
      return Vals[ind];
    }

    T &operator[](Index ind) {
      assert(ind >= 0 && ind < N);
      return Vals[ind];
    }

    T &back() {
      return at(N - 1);
    }
    const T &back() const {
      return at(N - 1);
    }

    //Return pointer to the first element
    T const *begin() const {
      return Vals;
    }
    //Return pointer to the first element
    T const *cbegin() const {
      return Vals;
    }
    T *begin() {
      return Vals;
    }

    //Return pointer to first point in memory beyond the last element
    T const *end() const {
      return Vals + N;
    }
    //Return pointer to first point in memory beyond the last element
    T const *cend() const {
      return Vals + N;
    }
    T *end() {
      return Vals + N;
    }


    //MUTATORS
    void push_back(const T &toPush);

    void pop_back() {
      if(N) Vals[--N].~T();
    }
    void remove(Index ind);
    void clear() {
      while(N) Vals[--N].~T();
    }

    void resize(Index new_N);
    void resize(Index new_N, const T &fill_val);
    void reserve(Index new_max);

    template <typename CompareType>
    void sort(const CompareType &comp);
    void sort(Array<Index> &ind_order);
    void sort();
    Array &append(const Array &new_tail);
    Array &append_unique(const Array &new_tail);

    void swap_elem(Index i, Index j) {
      std::swap(at(i), at(j));
    };

    Array &permute(const Array<Index> &perm_array);
    Array &ipermute(const Array<Index> &perm_array);
    bool next_permute();

    ReturnArray<Index> as_perm_inverse() const;
    ReturnArray<Index> as_perm_transform_by(const Array<Index> &trans_perm) const;
    // INSPECTORS

    const T &max() const;
    const T &min() const;

    // copy portion of *this Array into a new Array
    ReturnArray<T> sub_array(Index ind_begin, Index ind_end)const;

    T sum() const;

    bool is_ascending() const;
    bool is_descending() const;
    bool is_constant() const;
    bool is_permute() const;
    bool has_fixed_points() const;

    // COMPARISONS

    bool operator==(const Array<T> &RHS) const;
    bool operator!=(const Array<T> &RHS) const {
      return !((*this) == RHS);
    }
    bool operator<(const Array<T> &RHS) const;
    bool operator>(const Array<T> &RHS) const;
    bool operator<=(const Array<T> &RHS) const {
      return !((*this) > RHS);
    }
    bool operator>=(const Array<T> &RHS) const {
      return !((*this) < RHS);
    }

    bool all_in(const Array &superset) const;
    Index coincidence(const Array &superset) const;
    Index incidences(const T &test_elem) const;
    Index find(const T &test_elem) const;
    ///Same as find, but starts from the last element of the Array
    Index reverse_find(const T &test_elem) const;
    Index almost_find(const T &test_elem, double tol_val = TOL) const;
    ///Same as almost_find, but start from the last element of the Array
    Index almost_reverse_find(const T &test_elem, double tol_val = TOL) const;
    bool contains(const T &test_elem) const {
      return find(test_elem) < N;
    }
    bool almost_contains(const T &test_elem, double tol_val = TOL) const {
      return almost_find(test_elem, tol_val) < N;
    }

    // I/O
    void print_column(std::ostream &stream, const std::string &indent = "")const;
  };

  template<typename T>
  jsonParser &to_json(const Array<T> &value, jsonParser &json) {
    json.put_array();
    for(Index i = 0; i < value.size(); i++)
      json.push_back(value[i]);
    return json;
  }

  /// This requires that 'T::T()' exists, if not, you must do this by hand
  template<typename T>
  void from_json(Array<T> &value, const jsonParser &json) {
    try {
      value.resize(json.size());
      for(int i = 0; i < json.size(); i++)
        from_json(value[i], json[i]);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename T>
  class ReturnArray : public Array<T> {
  public:
    using Array<T>::swap;

    ReturnArray(Array<T> &init_tens) {
      swap(init_tens);
    }

    ReturnArray &operator=(Array<T> &RHS) {
      swap(RHS);
      return *this;
    }

  };

  //*******************************************************************************************

  template<typename T>
  Array<T>::Array(ReturnArray<T> &RHS) : N(0), NMax(0), Vals(nullptr) {
    swap(RHS);
  }

  //*******************************************************************************************
  /// Returns an array with the sequence (initial, ++initial, ..., final), inclusive
  template<typename T>
  ReturnArray<T> Array<T>::sequence(const T &initial, const T &final) {
    Array<T> seq;
    if(!(initial <= final))
      return seq;

    seq.push_back(initial);

    while(seq.back() < final) {
      seq.push_back(seq.back());
      seq.back()++;
    }
    return seq;
  }

  //*******************************************************************************************
  /// Returns an array with the sequence (initial, initial+increment, ..., final?),
  /// inclusive if final is in the sequence
  template<typename T>
  ReturnArray<T> Array<T>::sequence(const T &initial, const T &increment, const T &final) {
    Array<T> seq;
    seq.push_back(initial);
    seq.push_back(initial += increment);

    if(seq[0] < seq[1]) { // increasing case
      if(final < seq[1]) {
        seq.pop_back();
        if(final < initial) {
          seq.clear();
        }
        return seq;
      }

      while(seq.back() < final) {
        seq.push_back(seq.back());
        seq.back() += increment;
      }
      if(final < seq.back())
        seq.pop_back();
      return seq;
    }
    else if(seq[1] < seq[0]) { // decreasing case
      if(seq[1] < final) {
        seq.pop_back();
        if(initial < final) {
          seq.clear();
        }
        return seq;
      }

      while(final < seq.back()) {
        seq.push_back(seq.back());
        seq.back() += increment;
      }
      if(seq.back() < final)
        seq.pop_back();
      return seq;
    }
    else { // increment=0 case
      seq.pop_back();
      if(initial < final || final < initial)
        seq.clear();
      return seq;
    }

  }

  //*******************************************************************************************

  template<typename T>
  Array<T> &Array<T>::operator=(const Array<T> &RHS) {
    if(this == &RHS) {
      return *this;
    }

    Index i;
    if(NMax < RHS.size()) {
      clear();
      reserve(RHS.size());
      for(i = 0; i < RHS.size(); i++)
        push_back(RHS[i]);
      return *this;
    }

    while(N > RHS.size())
      pop_back();

    i = N;

    while(N < RHS.size())
      push_back(RHS[N]);

    while(i--)
      Vals[i] = RHS[i];

    return *this;
  }

  //*******************************************************************************************

  template<typename T>
  Array<T> &Array<T>::operator=(ReturnArray<T> &RHS) {
    swap(RHS);
  }

  //*******************************************************************************************
  template<typename T>
  void Array<T>::swap(Array<T> &RHS) {
    if(this == &RHS) return;

    Index tN(N), tNMax(NMax);
    T *tVals(Vals);

    Vals = RHS.Vals;
    N = RHS.N;
    NMax = RHS.NMax;

    RHS.Vals = tVals;
    RHS.N = tN;
    RHS.NMax = tNMax;

    return;
  }

  //*******************************************************************************************
  template<typename T>
  void Array<T>::resize(Index new_N) {
    clear();
    if(new_N > 0) {
      reserve(new_N);
      while(N < new_N)
        new(Vals + (N++)) T();
    }
    return;
  }

  //*******************************************************************************************

  template<typename T>
  void Array<T>::resize(Index new_N, const T &fill_val) {
    clear();
    reserve(new_N);
    while(N < new_N)
      new(Vals + (N++)) T(fill_val);
    return;
  }

  //*******************************************************************************************
  template<typename T>
  void Array<T>::reserve(Index new_max) {

    if(new_max <= N) return;

    T *tVal(nullptr);
    if(new_max) {
      tVal = static_cast<T *>(operator new(new_max * sizeof(T)));
      for(Index i = 0; i < N; i++) {
        new(tVal + i) T(Vals[i]);
        Vals[i].~T();
      }
    }
    if(Vals)
      operator delete(Vals);
    Vals = tVal;
    NMax = new_max;
    return;
  }

  //*******************************************************************************************

  template<typename T>
  void Array<T>::push_back(const T &toPush) {
    //If N==NMax, we must reallocate memory.  This is not done via Array::reserve.
    //Instead, we do it within Array::push_back
    if(N == NMax) {
      Index new_Max = (N * ARRAY_EXTENSION_FACTOR() > N + ARRAY_MIN_EXTRA_SPACE()) ?
                      (Index)(N * ARRAY_EXTENSION_FACTOR()) : N + ARRAY_MIN_EXTRA_SPACE();

      T *tVal(nullptr);

      tVal = static_cast<T *>(operator new(new_Max * sizeof(T)));

      //first, add the new element; this prevents aliasing problems
      new(tVal + N) T(toPush);
      for(Index i = 0; i < N; i++) {
        new(tVal + i) T(Vals[i]);
        Vals[i].~T();
      }

      if(Vals)
        operator delete(Vals);
      Vals = tVal;
      NMax = new_Max;

      N++;
    }
    else {
      new(Vals + (N++)) T(toPush);
    }

  }

  //*******************************************************************************************
  template<typename T>
  void Array<T>::remove(Index ind) {

    for(Index i = ind + 1; i < N; i++)
      at(i - 1) = at(i);

    Vals[--N].~T();

    return;
  }

  //*******************************************************************************************
  template<typename T>
  bool Array<T>::operator==(const Array<T> &RHS) const {
    if(N != RHS.N)
      return false;

    Index i;
    for(i = 0; i < N; i++) {
      if(!(at(i) == RHS[i])) {
        return false;
      }
    }

    return true;
  }

  //*******************************************************************************************
  template<typename T>
  bool Array<T>::operator<(const Array<T> &RHS) const {
    if(size() != RHS.size())
      return size() < RHS.size();

    for(Index i = 0; i < size(); i++) {
      if(at(i) < RHS[i])
        return true;
      else if(at(i) > RHS[i])
        return false;
    }
    return false;
  }

  //*******************************************************************************************
  template<typename T>
  bool Array<T>::operator>(const Array<T> &RHS) const {
    if(size() != RHS.size())
      return size() > RHS.size();

    for(Index i = 0; i < size(); i++) {
      if(at(i) > RHS[i])
        return true;
      else if(at(i) < RHS[i])
        return false;
    }
    return false;
  }

  //*******************************************************************************************
  template<typename T>
  const T &Array<T>::max() const {
    if(!size()) {
      std::cerr << "ERROR: Tried to find maximum value of empty array! Exiting...\n";
      assert(0);
      exit(1);
    }

    Index imax = 0;

    for(Index i = 1; i < size(); i++) {
      if(at(imax) < at(i))
        imax = i;
    }
    return at(imax);
  }

  //*******************************************************************************************

  template<typename T>
  const T &Array<T>::min() const {
    if(!size()) {
      std::cerr << "ERROR: Tried to find minimum value of empty array! Exiting...\n";
      assert(0);
      exit(1);
    }

    Index imin = 0;
    for(Index i = 1; i < size(); i++) {
      if(at(i) < at(imin))
        imin = i;
    }
    return at(imin);
  }

  //*******************************************************************************************

  template<typename T>
  T Array<T>::sum() const {
    T result(0);
    for(Index i = 0; i < size(); i++) {
      result += at(i);

    }
    return result;
  }

  //*******************************************************************************************
  template<typename T>
  ReturnArray<T> Array<T>::sub_array(Index ind_begin, Index ind_end)const {
    assert(ind_begin <= ind_end && ind_end < size() && "Array::operator() indices out of bounds!");
    Array<T> sub;
    sub.reserve(ind_end - ind_begin);
    while(ind_begin <= ind_end) {
      sub.push_back(at(ind_begin++));
    }
    return sub;
  }
  //*******************************************************************************************

  template<typename T>
  bool Array<T>::all_in(const Array &superset) const {
    for(Index i = 0; i < N; i++) {
      if(!superset.contains(at(i))) {
        return false;
      }
    }

    return true;
  }

  //*******************************************************************************************
  //returns number of elements of *this that are coincident with 'superset',
  // neglecting order.

  template<typename T>
  Index Array<T>::coincidence(const Array &superset) const {
    Index nco(0);
    for(Index i = 0; i < N; i++) {
      if(superset.contains(at(i))) {
        nco++;
      }
    }

    return nco;
  }

  //*******************************************************************************************

  template< typename T>
  Index Array<T>::incidences(const T &test_elem) const {
    Index ni(0);
    for(Index i = 0; i < N; i++) {
      if(at(i) == test_elem) {
        ni++;
      }
    }

    return ni;
  }

  //*******************************************************************************************

  template< typename T>
  Index Array<T>::find(const T &test_elem) const {
    for(Index i = 0; i < N; i++) {
      if(at(i) == test_elem) {
        return i;
      }
    }

    return N;
  }

  //*******************************************************************************************

  template< typename T>
  Index Array<T>::reverse_find(const T &test_elem) const {
    for(Index i = N - 1; i >= 0; i--) {
      if(at(i) == test_elem) {
        return i;
      };
    };

    return N;
  }

  //************************************************************

  template< typename T>
  Index Array<T>::almost_find(const T &test_elem, double tol_val) const {
    for(Index i = 0; i < N; i++) {
      if(almost_equal(at(i), test_elem, tol_val)) {
        return i;
      }
    }

    return N;
  }

  //*******************************************************************************************
  template< typename T>
  void Array<T>::print_column(std::ostream &stream, const std::string &indent)const {
    for(Index i = 0; i < N; i++)
      stream << indent << at(i) << std::endl;
  }
  //*******************************************************************************************

  template< typename T>
  Index Array<T>::almost_reverse_find(const T &test_elem, double tol_val) const {
    for(Index i = N - 1; i >= 0; i--) {
      if(almost_equal(at(i), test_elem, tol_val)) {
        return i;
      };
    };

    return N;
  }

  //************************************************************

  template< typename T>
  template <typename CompareType>
  void Array<T>::sort(const CompareType &comp) {
    /// quicksort sorting algorithm
    ///  - assumes that CompareType::compare(T a_thing, T b_thing) exists
    ///  - End results is that CompareType::compare(at(i), at(j)) is true for all i<j
    // Adapted from implementation by Brian Puchala

    int left, right, middle, pivot, curr, i;
    Array<int> queue_left, queue_right;

    // estimate queue size
    int cap = (int) 3 * (log(size()) / log(2));
    //cout << "cap: " << cap << endl;
    queue_left.reserve(cap);
    queue_right.reserve(cap);

    queue_left.push_back(0);
    queue_right.push_back(size() - 1);

    Index max_q = 0;


    do {
      left = queue_left[0];
      right = queue_right[0];
      queue_left.swap_elem(0, queue_left.size() - 1);
      queue_left.pop_back();

      queue_right.swap_elem(0, queue_right.size() - 1);
      queue_left.pop_back();

      if(left < right) {

        // avoid using (right+left) since this could be larger than max allowed value
        middle = left + (right - left) / 2;

        // choose median of (left, middle, right) for pivot
        if(comp.compare(at(left), at(middle))) {
          if(comp.compare(at(middle), at(right))) pivot = middle;
          else {
            if(comp.compare(at(left), at(right))) pivot = right;
            else pivot = left;
          }
        }
        else {
          if(comp.compare(at(right), at(middle))) pivot = middle;
          else {
            if(comp.compare(at(right), at(left))) pivot = right;
            else pivot = left;
          }
        }


        // move pivot to end
        swap_elem(pivot, right);

        // in the current region (left:right-1), seperate the values less than the pivot value (which is now at 'right')
        curr = left;
        for(i = left; i < right; i++) {
          if(comp.compare(at(i) , at(right))) {
            swap_elem(curr, i);
            curr++;
          }
        }

        // move the pivot value to the 'curr' position
        swap_elem(curr, right);

        // now everything in 'left:curr-1' is < *val[curr], and everything in 'curr+1:right' is >= *val[curr]
        //   so add those two regions to the queue

        if(curr != 0) {
          if(left != curr - 1) {
            //cout << "  add region: " << left << ":" << curr-1 << endl;
            queue_left.push_back(left);
            queue_right.push_back(curr - 1);
          }
        }

        if(curr + 1 != right) {
          //cout << "  add region: " << curr+1 << ":" << right << endl;
          queue_left.push_back(curr + 1);
          queue_right.push_back(right);
        }
        // repeat until there is nothing left
      }

      if(queue_left.size() > max_q)
        max_q = queue_left.size();


    }
    while(queue_left.size() > 0);

    //cout << "max_q: " << max_q << endl;
    //cout << "max_q/cap: " << (1.0*max_q)/(1.0*cap) << endl;
  }

  //*******************************************************************************************

  template<typename T>
  void Array<T>::sort(Array<Index> &ind_order) {
    ind_order.clear();
    for(Index i = 0; i < size(); i++)
      ind_order.push_back(i);

    for(Index i = 0; i < size(); i++) {
      for(Index j = i + 1; j < size(); j++) {
        if(at(j) < at(i)) {
          swap_elem(i, j);
          ind_order.swap_elem(i, j);
        }
      }
    }

  }

  //*******************************************************************************************

  template<typename T>
  void Array<T>::sort() {
    for(Index i = 0; i < size(); i++) {
      for(Index j = i + 1; j < size(); j++) {
        if(at(j) < at(i)) {
          swap_elem(i, j);
        }
      }
    }

  }
  //*******************************************************************************************
  template<typename T>
  Array<T> &Array<T>::append(const Array<T> &new_tail) {
    reserve(size() + new_tail.size());
    Index Nadd(new_tail.size());
    for(Index i = 0; i < Nadd; i++)
      push_back(new_tail[i]);
    return *this;
  }

  //*******************************************************************************************

  template<typename T>
  Array<T> &Array<T>::append_unique(const Array<T> &new_tail) {
    Index Nadd(new_tail.size());
    for(Index i = 0; i < Nadd; i++) {
      if(!contains(new_tail[i]))
        push_back(new_tail[i]);
    }

    return *this;
  }

  //*******************************************************************************************
  //perm_array gives order of indices after permutation in terms of original indices.
  //e.g., for my_array={2,3,5,8} and perm_array={3,2,1,4}
  //my_array.permute(perm_array)={5,3,2,8}
  template<typename T>
  Array<T> &Array<T>::permute(const Array<Index> &perm_array) {
    assert(perm_array.size() == size());
    Array<T> after_array;
    after_array.reserve(size());
    for(Index i = 0; i < perm_array.size(); i++)
      after_array.push_back(at(perm_array[i]));

    swap(after_array);
    return *this;
  }

  //*******************************************************************************************
  // generate inversely permuted copy of *this
  // perm_array gives order of indices after permutation in terms of original indices.
  // e.g., for my_array={2,3,5,8} and perm_array={3,2,1,4}
  // my_array.permute(perm_array)={5,3,2,8}
  template<typename T>
  Array<T> &Array<T>::ipermute(const Array<Index> &perm_array) {
    assert(perm_array.size() == size());
    Array<T> after_array(*this);

    for(Index i = 0; i < perm_array.size(); i++) {
      if(i != perm_array[i]) {
        after_array[perm_array[i]] = at(i);
      }
    }
    swap(after_array);
    return *this;
  }

  //*******************************************************************************************

  template<typename T>
  bool Array<T>::next_permute() {
    Index i(N - 2), j(N - 1);
    bool is_valid = true;
    while(valid_index(i) && at(i + 1) <= at(i))
      i--;

    if(valid_index(i)) {
      while(at(j) <= at(i))
        j--;

      swap_elem(i, j);
      i++;
      j = N - 1;
    }
    else {
      is_valid = false;
      i = 0;
    }

    while(i < j)
      swap_elem(i++, j--);

    return is_valid;
  }

  //*******************************************************************************************

  /// Construct permutation that undoes the permutation performed by 'this'
  /// Inverse operation is calculated the same, regardless of permutation convention
  template<typename T>
  ReturnArray<Index> Array<T>::as_perm_inverse() const {
    ReturnArray<Index> iperm_array(size(), 0);
    for(Index i = 0; i < size(); i++) {
      iperm_array[at(i)] = i;
    }
    return iperm_array;
  }

  //*******************************************************************************************
  /// Given N distinct objects labeled from 0 to N-1,
  /// a permutation 'P_permute' that physically permutes the objects (with labels) in terms of their labels,
  /// and a permutation 'L_permute' that permutes their labels only,
  /// rewrite 'P_permute' in terms of the relabeling induced by 'L_permute'
  /// Rearrange 'this' permutation to form an equivalent permutation for
  /// any list that has already been permuted by trans_perm.
  /// Does not (nearly certain of this) depend on permutation convention

  template<typename T>
  ReturnArray<Index> Array<T>::as_perm_transform_by(const Array<Index> &trans_perm) const {
    return ((trans_perm.as_perm_inverse()).permute(*this)).permute(trans_perm);
  }

  //*******************************************************************************************

  template<typename T>
  bool Array<T>::is_ascending() const {
    for(Index i = 0; i < N - 1; i++) {
      if(at(i + 1) < at(i))
        return false;
    }

    return true;
  }

  //*******************************************************************************************

  template<typename T>
  bool Array<T>::is_descending() const {
    for(Index i = 0; i < N - 1; i++) {
      if(at(i) < at(i + 1))
        return false;
    }

    return true;
  }

  //*******************************************************************************************

  template< typename T>
  bool Array<T>::is_constant() const {
    for(Index i = 0; i < N - 1; i++) {
      if(!(at(i) == at(i + 1)))
        return false;
    }

    return true;
  }

  //*******************************************************************************************

  /// Checks that Array contains values from 0 to perm_array.size()-1 and that no value is repeated
  /// does not depend on definition of permutation convention
  template< typename T>
  bool Array<T>::is_permute() const {
    for(Index i = 0; i < size(); i++) {
      if(at(i) < 0 || at(i) >= size() || find(at(i)) < i)
        return false;
    }
    return true;
  }

  //*******************************************************************************************

  /// Checks whether any values are equal to their index -- only valid for Array<Index>
  template< typename T>
  bool Array<T>::has_fixed_points() const {
    for(Index i = 0; i < size(); i++) {
      if(at(i) == i)
        return true;
    }
    return false;
  }


  //*******************************************************************************************
  template<class T>
  Array<T> array_cat(const Array<T> &A1, const Array<T> &A2) {
    return Array<T>(A1).append(A2);

  }

  //*******************************************************************************************
  template<class T>
  std::ostream &operator<<(std::ostream &out, const Array<T> &array_out) {
    if(array_out.size() == 0)
      out << "[empty]  ";
    for(Index i = 0; i < array_out.size(); i++) {
      out << array_out[i] << "  ";
    }
    return out;
  }



  //*******************************************************************************************
  /*
  template<class T>
  std::ostream &operator<<(std::ostream &out, const Array<T *> &array_out) {
    if(array_out.size() == 0)
      out << "[empty]  ";
    for(Index i = 0; i < array_out.size(); i++) {
      out << *array_out[i] << "  ";
    }
    return out;
  }
  */

  //*******************************************************************************************

  template<typename T>
  void swap(Array<T> &A1, Array<T> &A2) {
    A1.swap(A2);
    return;
  }

  //*******************************************************************************************

  class ArraySizeLessThan {
  public:
    template < typename T>
    bool compare(const Array<T> &a, const Array<T> &b) const {
      return a.size() < b.size();
    }
  };

  /** @} */

}

#endif /* ARRAY_H_ */
