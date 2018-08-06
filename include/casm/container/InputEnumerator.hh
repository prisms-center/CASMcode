#ifndef CASM_InputEnumerator
#define CASM_InputEnumerator

#include <string>

#include "casm/container/Enumerator.hh"

namespace CASM {

  template<typename ValueType, bool IsConst>
  class InputEnumeratorBase;

  template<typename ValueType, bool IsConst>
  class InputEnumIterator;


  template<typename ValueType, bool IsConst = true>
  class InputEnumIteratorBase : public ValEnumIterator<ValueType, IsConst> {

    friend InputEnumIterator<ValueType, IsConst>;

  public:

    using typename ValEnumIterator<ValueType, IsConst>::step_type;
    using typename ValEnumIterator<ValueType, IsConst>::value_type;
    using typename ValEnumIterator<ValueType, IsConst>::reference;


    InputEnumIteratorBase() {}

    InputEnumIteratorBase(InputEnumeratorBase<ValueType, IsConst> &enumerator, bool is_end) :
      ValEnumIterator<ValueType, IsConst>(enumerator), m_constructed_as_end(is_end) {}

    virtual ~InputEnumIteratorBase() {}


    using EnumIteratorBase::source;
    using EnumIteratorBase::name;

    /// Return current step number
    ///
    /// - Only valid if iterator refers to valid object (not end)
    virtual step_type step() const override {
      this->_assert_ptr();
      return _enum().step();
    }

    /// Check if end iterator
    ///
    /// Returns true if:
    /// - Iterator was constructed as end iterator
    /// - Or, !InputEnumeratorBase->valid()
    ///   - This means that the derived InputEnumerator class implementation
    ///     must call _invalidate in the increment implementation if the end
    ///     has been reached
    virtual bool is_end() const override {
      this->_assert_ptr();
      return _constructed_as_end() || !this->_enum().valid();
    }

    std::unique_ptr<InputEnumIteratorBase> clone() const {
      return std::unique_ptr<InputEnumIteratorBase>(this->_clone());
    }

  protected:
    using ValEnumIterator<ValueType, IsConst>::_enum_ptr;

    bool _constructed_as_end() const {
      return m_constructed_as_end;
    }

  private:

    virtual InputEnumIteratorBase<ValueType, IsConst> *_clone() const override {
      return new InputEnumIteratorBase<ValueType, IsConst>(*this);
    }

    /// boost::iterator_facade implementation
    virtual reference dereference() const override {
      return _enum().current();
    }

    /// boost::iterator_facade implementation
    virtual void increment() {
      _enum().increment();
    }

    /// static_cast enumerator pointer to InputEnumeratorBase&
    InputEnumeratorBase<ValueType, IsConst> &_enum() {
      return static_cast<InputEnumeratorBase<ValueType, IsConst>&>(*_enum_ptr());
    }

    /// static_cast enumerator pointer to InputEnumeratorBase&
    InputEnumeratorBase<ValueType, IsConst> &_enum() const {
      return static_cast<InputEnumeratorBase<ValueType, IsConst>&>(*_enum_ptr());
    }

    bool m_constructed_as_end;
  };


  template<typename ValueType, bool IsConst = true>
  class InputEnumIterator :

    public boost::iterator_facade <
    InputEnumIterator<ValueType, IsConst>,
    ValueType,
    std::input_iterator_tag,
    CASM_TMP::ConstSwitch<IsConst, ValueType> &,
    long > {

  public:

    typedef typename InputEnumIteratorBase<ValueType, IsConst>::step_type step_type;
    typedef ValueType value_type;
    typedef typename InputEnumIteratorBase<ValueType, IsConst>::reference reference;


    /// Default constructor
    InputEnumIterator() {}

    /// Construct iterator
    InputEnumIterator(const InputEnumIteratorBase<ValueType, IsConst> &it) :
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


  private:

    friend boost::iterator_core_access;

    /// boost::iterator_facade implementation
    void increment() {
      m_ptr->increment();
    }

    /// boost::iterator_facade implementation
    reference dereference() const {
      return m_ptr->dereference();
    }

    /// boost::iterator_facade implementation
    bool equal(const InputEnumIterator &B) const {
      return m_ptr->equal(*(B.m_ptr));
    }

    notstd::cloneable_ptr<InputEnumIteratorBase<ValueType, IsConst> > m_ptr;
  };


  /// \brief Base class for implementing specialized input enumerators
  ///
  /// InputEnumerator allow single pass iteration over some objects of ValueType.
  /// The specialized class should document any guaranteed properties of the
  /// enumerated objects. For instance, some enumerators may guarantee that
  /// objects will be unique, others may not. Some may guarantee objects will be
  /// primitive or in canonical form, others may not. Etc.
  ///
  /// Derived classes (ex. MyDerivedEnumClass) must implement:
  /// - constructors
  ///   - Must call _initialize or _set_current_ptr to set pointer to first valid object
  ///   - Ensure that the current step value = 0
  ///   - If due to the constructor arguments there are no valid objects, then
  ///     call _invalidate
  /// - virtual destructor
  /// - void increment()
  ///   - Implementation must do work, using _set_current_ptr if necessary to
  ///     point the enumerator at the object to be deferenced
  ///   - Then if still valid call _increment_step,
  ///     or if no longer valid ('end' has been reached) call _invalidate
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
  class InputEnumeratorBase : public ValEnumerator<ValueType, IsConst> {

    friend InputEnumIteratorBase<ValueType, IsConst>;

  public:

    using typename ValEnumerator<ValueType, IsConst>::step_type;
    using typename ValEnumerator<ValueType, IsConst>::value_type;
    using typename ValEnumerator<ValueType, IsConst>::reference;
    typedef InputEnumIterator<ValueType, IsConst> iterator;


    InputEnumeratorBase() {}

    virtual ~InputEnumeratorBase() {}


    using ValEnumerator<ValueType, IsConst>::step;
    using ValEnumerator<ValueType, IsConst>::valid;
    using ValEnumerator<ValueType, IsConst>::source; // virtual
    using ValEnumerator<ValueType, IsConst>::name; // pure virtual
    using ValEnumerator<ValueType, IsConst>::current;

    iterator begin() {
      return InputEnumIteratorBase<ValueType, IsConst>(*this, false);
    }

    iterator end() {
      return InputEnumIteratorBase<ValueType, IsConst>(*this, true);
    }


  protected:

    using ValEnumerator<ValueType, IsConst>::_initialize;
    using ValEnumerator<ValueType, IsConst>::_set_step;
    using ValEnumerator<ValueType, IsConst>::_increment_step;
    using ValEnumerator<ValueType, IsConst>::_invalidate;
    using ValEnumerator<ValueType, IsConst>::_set_current_ptr;
    using ValEnumerator<ValueType, IsConst>::_current;


  private:

    //using ValEnumerator<ValueType, IsConst>::_decrement_step; // hide
    virtual void increment() = 0;

  };

  /// \brief Generic input enumerator
  ///
  /// - Wraps a InputEnumeratorBase derived enumerator, so it can provide
  ///   single pass iteration using either an input enumerator or a random access
  ///   enumerator
  template<typename ValueType, bool IsConst = true>
  class InputEnumerator {

  public:

    typedef typename InputEnumeratorBase<ValueType, IsConst>::step_type step_type;
    typedef typename InputEnumeratorBase<ValueType, IsConst>::value_type value_type;
    typedef typename InputEnumeratorBase<ValueType, IsConst>::reference reference;
    typedef typename InputEnumeratorBase<ValueType, IsConst>::iterator iterator;


    InputEnumerator() {}

    InputEnumerator(std::unique_ptr<InputEnumeratorBase<ValueType, IsConst> > e) :
      ptr(notstd::clone(e)) {}

    InputEnumerator(const InputEnumerator &) = delete;
    InputEnumerator &operator=(const InputEnumerator &) = delete;


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

    /// Access the current ObjectType by reference
    reference current() const {
      return ptr->current();
    }

    iterator begin() {
      return ptr->begin();
    }

    iterator end() {
      return ptr->end();
    }

    std::unique_ptr<InputEnumeratorBase<ValueType, IsConst> > ptr;
  };

}

#endif
