#ifndef CASM_RandomAccessEnumerator
#define CASM_RandomAccessEnumerator

#include <string>

#include "casm/container/Enumerator.hh"
#include "casm/container/InputEnumerator.hh"

namespace CASM {

  template<typename ValueType, bool IsConst>
  class RandomAccessEnumIterator;

  template<typename ValueType, bool IsConst>
  class RandomAccessEnumeratorBase;

  template<typename ValueType, bool IsConst>
  class RandomAccessEnumerator;


  /// \brief RandomAccessEnumIterator implemenation
  ///
  /// - Holds a member variable 'm_it_step' indicating current 'step' value.
  /// - If 'step' == enumerator 'size', then the iterator has reached 'end'
  /// - To dereference an object, ensures that the enumerator being pointed at is
  ///   at step 'm_it_step' by calling 'goto_step'. Then the enumerator's 'current'
  ///   object reference can be safely returned.
  /// - Inherits from InputEnumIteratorBase so that it can be held by InputEnumIterator
  ///
  template<typename ValueType, bool IsConst = true>
  class RandomAccessEnumIteratorBase : public InputEnumIteratorBase<ValueType, IsConst> {

    friend RandomAccessEnumIterator<ValueType, IsConst>;

  public:

    using typename InputEnumIteratorBase<ValueType, IsConst>::step_type;
    using typename InputEnumIteratorBase<ValueType, IsConst>::value_type;
    using typename InputEnumIteratorBase<ValueType, IsConst>::reference;
    typedef step_type difference_type;


    RandomAccessEnumIteratorBase() {}

    RandomAccessEnumIteratorBase(RandomAccessEnumeratorBase<ValueType, IsConst> &enumerator, bool is_end) :
      InputEnumIteratorBase<ValueType, IsConst>(enumerator, is_end) {

      if(!is_end) {
        m_it_step = 0;
      }
      else {
        m_it_step = _enum().size();
      }
    }


    using InputEnumIteratorBase<ValueType, IsConst>::source;
    using InputEnumIteratorBase<ValueType, IsConst>::name;

    /// Return current step number
    ///
    /// - Only valid if iterator refers to valid object (not end)
    step_type step() const override {
      return m_it_step;
    }

    /// Check if end iterator
    ///
    /// Returns true if:
    /// - Iterator was constructed as end iterator
    /// - Or, iterator step value equals enumerator size
    bool is_end() const override {
      this->_assert_ptr();
      return _constructed_as_end() || step() == _enum().size();
    }

    std::unique_ptr<RandomAccessEnumIteratorBase> clone() const {
      return std::unique_ptr<RandomAccessEnumIteratorBase>(this->_clone());
    }

  protected:

    using InputEnumIteratorBase<ValueType, IsConst>::_enum_ptr;
    using InputEnumIteratorBase<ValueType, IsConst>::_constructed_as_end;

  private:

    /// boost::iterator_facade implementation
    reference dereference() const override {
      _enum().goto_step(m_it_step);
      return _enum().current();
    }

    /// boost::iterator_facade implementation
    void increment() override {
      ++m_it_step;
    }

    /// boost::iterator_facade implementation
    void decrement() {
      --m_it_step;
    }

    /// boost::iterator_facade implementation
    void advance(step_type n) {
      m_it_step += n;
    }

    /// boost::iterator_facade implementation
    difference_type distance_to(const RandomAccessEnumIteratorBase &B) const {
      return B.m_it_step - m_it_step;
    }

    /// clone implementation
    RandomAccessEnumIteratorBase<ValueType, IsConst> *_clone() const override {
      return new RandomAccessEnumIteratorBase<ValueType, IsConst>(*this);
    }

    /// static_cast enumerator pointer to RandomAccessEnumeratorBase&
    RandomAccessEnumeratorBase<ValueType, IsConst> &_enum() {
      return static_cast<RandomAccessEnumeratorBase<ValueType, IsConst>&>(*_enum_ptr());
    }

    /// static_cast enumerator pointer to RandomAccessEnumeratorBase&
    RandomAccessEnumeratorBase<ValueType, IsConst> &_enum() const {
      return static_cast<RandomAccessEnumeratorBase<ValueType, IsConst>&>(*_enum_ptr());
    }

    step_type m_it_step;
  };


  /// \brief The iterator type for RandomAccessEnumerator
  ///
  template<typename ValueType, bool IsConst = true>
  class RandomAccessEnumIterator :

    public boost::iterator_facade <
    RandomAccessEnumIterator<ValueType, IsConst>,
    ValueType,
    boost::random_access_traversal_tag,
    CASM_TMP::ConstSwitch<IsConst, ValueType> &,
    long > {

  public:

    typedef typename RandomAccessEnumIteratorBase<ValueType, IsConst>::step_type step_type;
    typedef ValueType value_type;
    typedef typename RandomAccessEnumIteratorBase<ValueType, IsConst>::reference reference;
    typedef step_type difference_type;


    /// Default constructor
    RandomAccessEnumIterator() {}

    /// Construct iterator
    RandomAccessEnumIterator(const RandomAccessEnumIteratorBase<ValueType, IsConst> &it) :
      m_ptr(notstd::clone(it)) {}


    step_type step() const {
      return m_ptr->step();
    }

    jsonParser source() const {
      return m_ptr->source();
    }

    std::string name() const {
      return m_ptr->name();
    }

    std::unique_ptr<RandomAccessEnumIterator> clone() const {
      return std::unique_ptr<RandomAccessEnumIterator>(this->_clone());
    }


  private:

    friend boost::iterator_core_access;

    /// clone implementation
    RandomAccessEnumIterator *_clone() const {
      return new RandomAccessEnumIterator(*this);
    }

    /// boost::iterator_facade implementation
    bool equal(const RandomAccessEnumIterator &B) const {
      return m_ptr->equal(*(B.m_ptr));
    }

    reference dereference() const {
      return m_ptr->dereference();
    }

    void increment() {
      m_ptr->increment();
    }

    void decrement() {
      m_ptr->decrement();
    }

    void advance(step_type n) {
      m_ptr->advance(n);
    }

    difference_type distance_to(const RandomAccessEnumIterator &B) const {
      return m_ptr->distance_to(*(B.m_ptr));
    }

    notstd::cloneable_ptr<RandomAccessEnumIteratorBase<ValueType, IsConst> > m_ptr;
  };


  /// \brief Base class for implementing specialized random access enumerators
  ///
  /// InputEnumerator allow random access iteration over some objects of ValueType.
  /// The specialized class should document any guaranteed properties of the
  /// enumerated objects. For instance, some enumerators may guarantee that
  /// objects will be unique, others may not. Some may guarantee objects will be
  /// primitive or in canonical form, others may not. Etc.
  ///
  /// - Inherits from InputEnumeratorBase so that InputEnumerator can hold
  ///   both input enumerators and random access enumerators
  ///
  /// Derived classes (ex. MyDerivedEnumClass) must implement:
  /// - constructors
  ///   - Must call _initialize or _set_current_ptr to set pointer to first valid object
  ///   - Ensure that the current step value = 0, via _initialize or _set_step
  ///   - Ensure that the size of the enumerator is set, via _set_size
  ///   - If due to the constructor arguments there are no valid objects, then
  ///     call _invalidate
  /// - virtual destructor
  /// - void goto_step(step_type n)
  ///   - Should do work to determine the 'n'-th object and then call
  ///     goto_step(step_type n, value_type* ptr) setting a pointer to that object
  ///   - If 'n' equals the current step, can immediately return
  ///
  /// Derived classes should contain ENUMERATOR_MEMBERS(MyDerivedClass) in the
  /// class definition to automatically implement:
  ///   - name
  ///
  /// In CASM namespace near the derived class definition, include:
  /// - For enumerators only meant to be used internally:
  ///   - In header:
  ///     - ENUMERATOR_TRAITS(MyDerivedEnumClass)
  ///   - In the source code file:
  ///     - const std::string CASM_TMP::traits<MyDerivedEnumClass>::name = "MyDerivedEnumClass";
  /// - For enumerators only meant to be added to the API:
  ///   - In header:
  ///     - ENUMERATOR_INTERFACE_TRAITS(MyDerivedClass)
  ///   - In the source code file:
  ///     - \code
  ///       const std::string CASM_TMP::traits<MyDerivedEnumClass>::name = "MyDerivedEnumClass";
  ///       const std::string CASM_TMP::traits<MyDerivedEnumClass>::help =
  ///         "MyDerivedEnumClass: \n\n"
  ///
  ///         "  kwarg1: type (default=X) \n"
  ///         "    Description of kwarg1... \n\n"
  ///
  ///         "  kwarg2: type (default=X) \n"
  ///         "    Description of kwarg2... \n\n"
  ///
  ///         "  ... \n";
  ///
  ///       int EnumInterface<MyDerivedEnumClass>::run(PrimClex &primclex, const jsonParser &kwargs) const {
  ///         ...
  ///         implementation to run a MyDerivedEnumClass enumerator constructed
  ///         using the JSON input in 'kwargs' and store & save results in primclex
  ///         ...
  ///         return (returncode);
  ///       }
  ///       \endcode
  ///
  ///
  template<typename ValueType, bool IsConst = true>
  class RandomAccessEnumeratorBase : public InputEnumeratorBase<ValueType, IsConst> {

    friend RandomAccessEnumIteratorBase<ValueType, IsConst>;

  public:

    using typename InputEnumeratorBase<ValueType, IsConst>::step_type;
    using typename InputEnumeratorBase<ValueType, IsConst>::value_type;
    using typename InputEnumeratorBase<ValueType, IsConst>::reference;
    typedef RandomAccessEnumIterator<ValueType, IsConst> iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef step_type size_type;


    RandomAccessEnumeratorBase(size_type size = 0) :
      InputEnumeratorBase<ValueType, IsConst>(),
      m_size(size) {}

    ~RandomAccessEnumeratorBase() {}


    using InputEnumeratorBase<ValueType, IsConst>::step;
    using InputEnumeratorBase<ValueType, IsConst>::valid;
    using InputEnumeratorBase<ValueType, IsConst>::source; // virtual
    using InputEnumeratorBase<ValueType, IsConst>::name; // pure virtual
    using InputEnumeratorBase<ValueType, IsConst>::current;

    iterator begin() {
      return RandomAccessEnumIteratorBase<ValueType, IsConst>(*this, false);
    }

    iterator end() {
      return RandomAccessEnumIteratorBase<ValueType, IsConst>(*this, true);
    }

    reverse_iterator rbegin() {
      return reverse_iterator(end());
    }

    reverse_iterator rend() {
      return reverse_iterator(begin());
    }

    /// Number of elements in enumerator
    size_type size() const {
      return m_size;
    };

    /// Reference to first element in enumerator
    reference front() {
      this->goto_step(0);
      return this->dereference();
    }

    /// Reference to last element in enumerator
    reference back() {
      this->goto_step(size() - 1);
      return this->dereference();
    }

    /// Reference to first element in enumerator
    reference initial() {
      return front();
    }

    /// Reference to last element in enumerator
    reference final() {
      return back();
    }

    /// Reference an element in enumerator
    reference operator[](size_type n) {
      this->goto_step(n);
      return this->dereference();
    }

  protected:

    using ValEnumerator<ValueType, IsConst>::_initialize;
    using ValEnumerator<ValueType, IsConst>::_set_step;
    using ValEnumerator<ValueType, IsConst>::_increment_step;
    using ValEnumerator<ValueType, IsConst>::_decrement_step;
    using ValEnumerator<ValueType, IsConst>::_invalidate;
    using ValEnumerator<ValueType, IsConst>::_set_current_ptr;
    using ValEnumerator<ValueType, IsConst>::_current;

    /// Set size value
    void _set_size(size_type val) {
      m_size = val;
    }

  private:

    /// Must be implemented in derived enumerator
    ///
    /// - Should return a pointer to the object 'at' step n
    /// - May be implemented by returning pointers to different objects, or
    ///   by modifying the object being pointed at and returning the same pointer
    virtual value_type *at_step(step_type n) = 0;

    /// Sets step, current pointer via call to 'at_step', and enumerator validity
    /// - If n in valid range, calls at_step to set the current object and
    ///   sets the enumerator state to valid
    /// - If n not in valid range, will set enumerator state to not valid
    void goto_step(step_type n) {
      if(this->step() == n) {
        return;
      }
      this->_set_step(n);
      if(n >= 0 && n < this->size()) {
        this->_validate();
        this->_set_current_ptr(this->at_step(n));
      }
      else {
        this->_invalidate();
      }
    }

    /// Pure virtual function inherited from InputEnumeratorBase,
    ///   needed when incrementing from InputEnumeratorBase*
    void increment() override {
      goto_step(step() + 1);
    }

    size_type m_size;
  };


  /// \brief Generic random access enumerator
  ///
  /// - Wraps a RandomAccessEnumeratorBase derived enumerator
  template<typename ValueType, bool IsConst = true>
  class RandomAccessEnumerator {

  public:

    typedef typename RandomAccessEnumeratorBase<ValueType, IsConst>::step_type step_type;
    typedef ValueType value_type;
    typedef typename RandomAccessEnumeratorBase<ValueType, IsConst>::reference reference;
    typedef RandomAccessEnumIterator<ValueType, IsConst> iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef typename RandomAccessEnumeratorBase<ValueType, IsConst>::size_type size_type;

    RandomAccessEnumerator() {}

    RandomAccessEnumerator(std::unique_ptr<RandomAccessEnumeratorBase<ValueType, IsConst> > e) :
      ptr(e.release()) {}

    RandomAccessEnumerator(const RandomAccessEnumerator &) = delete;
    RandomAccessEnumerator &operator=(const RandomAccessEnumerator &) = delete;


    iterator begin() {
      return ptr->begin();
    }

    iterator end() {
      return ptr->end();
    }

    reverse_iterator rbegin() {
      return reverse_iterator(end());
    }

    reverse_iterator rend() {
      return reverse_iterator(begin());
    }

    /// Access the current ObjectType by reference
    reference current() const {
      return ptr->current();
    }

    /// Increments with each enumerated object
    step_type step() const {
      return ptr->step();
    }

    /// Returns false if enumeration is complete
    bool valid() const {
      return ptr->valid();
    }

    jsonParser source(step_type step) const {
      return ptr->source(step);
    }

    /// Derived enumerators must implement name
    std::string name() const {
      return ptr->name();
    }

    std::unique_ptr<RandomAccessEnumeratorBase<ValueType, IsConst> > ptr;
  };

}

#endif
