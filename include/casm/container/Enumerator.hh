#ifndef CASM_Enumerator
#define CASM_Enumerator

#include <string>

#include "casm/casm_io/jsonParser.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/misc/unique_cloneable_map.hh"

namespace CASM {

  /** \defgroup Enumerator
   *
   *  \ingroup Container
   *  \brief Algorithms for enumerating objects such as Lattice, Supercell, Configuration
   *
   *  Enumerators are classes that implement methods to enumerate Supercells,
   *  Configurations, and other objects by providing iterators over a range of
   *  enumerated objects. The objects are not usually stored in the enumerators
   *  but constructed as iterators are modified or dereferenced. New enumerators
   *  can be created by inheriting from:
   *  - InputEnumeratorBase (for single-pass  enumerators)
   *  - RandomAccessEnumeratorBase (for multi-pass, random-access enumerators)
   *
   *  Some example enumerators are:
   *  - ScelEnumT, ScelEnumByNameT, ScelEnumByPropsT
   *  - ConfigEnumAllOccupations
   *  - ConfigEnumInterpolation
   *  - ConfigEnumStrain
   *  - SuperConfigEnum
   *
   *  Enumerators are required to "know" their own name by implementing a
   *  traits class with 'name' as a const std::string member. For enumerators
   *  meant to be accessible via the 'casm enum' API, the traits class must also
   *  have a 'help' const std::string member explaining the enumerator and the
   *  input parameters, and a 'CASM::EnumInterface<EnumMethod>::run' method which
   *  implements collecting input parameters from the CLI, executing the enumerator
   *  and saving the resulting objects. Some helper functions exist to make
   *  parsing input and saving results easier:
   *
   *  - To collect CLI input and construct an ScelEnum which enumerates canonical
   *  Supercell: ::make_enumerator_scel_enum
   *  - To collect CLI input and construct an SupercellEnumerator<Lattice> which
   *  enumerates super-lattices which may not be canonical: ::make_enumerator_superlat_enum
   *  - To save results from enumerators of unique, primitive, canonical
   *  Configurations: ::insert_unique_canon_configs
   *  - To save results from enumerators of general Configurations: ::insert_configs
   *
   *
   *  For Enumerators only meant to be used internally the required members are:
   *
   *  Variables:
   *  - \code public: static const std::string enumerator_name; \endcode
   *
   *  Functions:
   *  - \code public: std::string name() const override { return enumerator_name; } \endcode
   *
   *
   *  For Enumerators to be added to the API the required members are:
   *
   *  Variables:
   *  - \code public: static const std::string enumerator_name; \endcode
   *  - \code public: static const std::string interface_help; \endcode
   *
   *  Functions:
   *  - public: std::string name() const override { return enumerator_name; } \endcode
   *  - public: static int run(PrimClex &primclex, const jsonParser &kwargs, const Completer::EnumOption &enum_opt); \endcode
   *
   *  To enable use as a plugin:
   *  - \code
   *    extern "C" {
   *      CASM::EnumInterfaceBase *make_EnumMethod_interface() {
   *        return new CASM::EnumInterface<CASM::ConfigEnumAllOccupations>();
   *      }
   *    }
   *    \endcode
   *  - To use an enumerator as a plugin for an existing CASM project, place the
   *    source code in the `.casm/enumerators` directory in a file named
   *    `EnumMethodName.cc`, where `EnumMethod` is the name of the enumerator
   *    class.

      @{
  */

  template<typename ValueType, bool IsConst>
  class ValEnumerator;

  // ---- Enumerator ---------------------

  /// \brief Abstract base class for enumerators
  ///
  /// - To implement a new enumeration method do not inherit from this
  ///   directly. Instead, inherit from either InputEnum (for single-pass
  ///   enumerators) or RandomAccessEnum (for multi-pass, random-access
  ///   enumerators)
  /// - The Enumerator base class only holds a pointer to the 'current'
  ///   object, the current 'step' index, and a 'valid' flag
  ///
  class EnumeratorBase {

  public:

    typedef long step_type;

    /// Default constructor
    EnumeratorBase() :
      m_valid(false),
      m_step(0) {}

    ~EnumeratorBase() {}

    /// Increments with each enumerated object
    step_type step() const {
      return m_step;
    }

    /// Returns false if enumeration is complete
    bool valid() const {
      return m_valid;
    }

    /// Default Object source just uses step#
    ///
    /// Returns:
    /// \code
    /// {
    ///   "enumerated_by": "<enumerator_type>",
    ///   "step": <step #>
    /// }
    /// \endcode
    virtual jsonParser source(step_type step) const {
      jsonParser src;
      src["enumerated_by"] = this->name();
      src["step"] = step;
      return src;
    }

    /// \brief Derived enumerators must implement name, via ENUM_MEMBERS
    virtual std::string name() const = 0;


  protected:

    /// Initialize
    ///
    /// - Sets step to 0
    /// - Sets valid to true
    void _initialize() {
      m_step = 0;
      m_valid = true;
    }

    /// Set current step value
    void _set_step(step_type val) {
      m_step = val;
    }

    /// Increment current step value
    void _increment_step() {
      ++m_step;
    }

    /// Decrement current step value
    void _decrement_step() {
      --m_step;
    }

    /// Call if enumeration complete
    void _invalidate() {
      m_valid = false;
    }

    /// Used if random access enumerator step is moved into valid range
    void _validate() {
      m_valid = true;
    }

  private:

    bool m_valid;

    step_type m_step;

  };

  template<typename ValueType, bool IsConst = true>
  class ValEnumerator : public EnumeratorBase {

  public:

    typedef ValueType value_type;
    typedef CASM_TMP::ConstSwitch<IsConst, ValueType> &reference;
    using EnumeratorBase::step_type;

    ValEnumerator():
      m_current_ptr(nullptr) {}

    virtual ~ValEnumerator() {}

    // -- from EnumeratorBase --

  public:

    using EnumeratorBase::step;
    using EnumeratorBase::valid;
    using EnumeratorBase::source; // virtual
    using EnumeratorBase::name; // pure virtual

  protected:

    /// Initialize
    ///
    /// - Sets current to point at _initial
    /// - Sets step to 0
    /// - Sets valid to true
    void _initialize(value_type *_initial) {
      _set_current_ptr(_initial);
      EnumeratorBase::_initialize();
    }

    using EnumeratorBase::_set_step;
    using EnumeratorBase::_increment_step;
    using EnumeratorBase::_decrement_step;
    using EnumeratorBase::_invalidate;


    // -- For ValEnumerator --

  public:

    /// Access the current ObjectType by reference
    reference current() const {
      return *m_current_ptr;
    }

  protected:

    /// Change the pointer
    void _set_current_ptr(value_type *_new) {
      m_current_ptr = _new;
    }

    /// Access the current ObjectType by reference
    value_type &_current() {
      return *m_current_ptr;
    }

  private:

    value_type *m_current_ptr;

  };


  class EnumIteratorBase {

  public:

    typedef EnumeratorBase::step_type step_type;

    /// Default Constructor
    EnumIteratorBase() :
      m_enum_ptr(nullptr) {}

    /// Constructor
    EnumIteratorBase(EnumeratorBase &enumerator) :
      m_enum_ptr(&enumerator) {}

    virtual ~EnumIteratorBase() {}

    /// Return current step number
    ///
    /// - Only valid if iterator refers to valid object (not end)
    virtual step_type step() const = 0;

    /// Uses 'step' and enumerator class 'source' implementation
    ///
    /// - Only valid if iterator refers to valid object (not end)
    jsonParser source() const {
      return _enum().source(this->step());
    }

    /// Uses enumerator class 'name' implementation
    std::string name() const {
      return m_enum_ptr->name();
    }

    /// Returns true if 'end' iterator
    virtual bool is_end() const = 0;

    std::unique_ptr<EnumIteratorBase> clone() const {
      return std::unique_ptr<EnumIteratorBase>(this->_clone());
    }


  protected:

    void _assert_same_ptr(const EnumIteratorBase &other) const {
      assert((m_enum_ptr == other.m_enum_ptr) &&
             "Comparing EnumIterator referring to different enumerators is not allowed!");
    }

    void _assert_ptr() const {
      assert(m_enum_ptr && "EnumIterator does not point to any enumerator!");
    }

    void _assert_valid() const {
      assert(m_enum_ptr->valid() && "EnumIterator points to an invalid enumerator!");
    }

    /// \brief boost::iterator_facade implementation
    ///
    /// - Uses 'is_end' implementation to check if both iterators are 'end'
    /// - If both are not end, then compares iterator 'step'
    bool equal(const EnumIteratorBase &other) const {
      bool this_is_end = this->is_end();
      bool other_is_end = other.is_end();

      if(this_is_end != other_is_end) {
        return false;
      }

      if(this_is_end) {
        return true;
      }

      return (this->step() == other.step());
    }

    EnumeratorBase *_enum_ptr() {
      return m_enum_ptr;
    }

    EnumeratorBase *_enum_ptr() const {
      return m_enum_ptr;
    }

  private:

    virtual EnumIteratorBase *_clone() const = 0;

    EnumeratorBase &_enum() {
      return *m_enum_ptr;
    }

    EnumeratorBase &_enum() const {
      return *m_enum_ptr;
    }

    // pointer to Enumerator
    EnumeratorBase *m_enum_ptr;

  };

  template<typename ValueType, bool IsConst = true>
  class ValEnumIterator : public EnumIteratorBase {

  public:

    using EnumIteratorBase::step_type;
    typedef ValueType value_type;
    typedef typename ValEnumerator<ValueType, IsConst>::reference reference;


    ValEnumIterator() {}

    ValEnumIterator(ValEnumerator<ValueType, IsConst> &enumerator) :
      EnumIteratorBase(enumerator) {}

    virtual ~ValEnumIterator() {}

    using EnumIteratorBase::step; // pure virtual
    using EnumIteratorBase::source;
    using EnumIteratorBase::name;
    using EnumIteratorBase::is_end; // pure virtual

    std::unique_ptr<EnumIteratorBase> clone() const {
      return std::unique_ptr<EnumIteratorBase>(this->_clone());
    }

  protected:

    using EnumIteratorBase::_assert_same_ptr;
    using EnumIteratorBase::_assert_ptr;
    using EnumIteratorBase::_assert_valid;
    using EnumIteratorBase::equal;
    using EnumIteratorBase::_enum_ptr;

  private:

    virtual EnumIteratorBase *_clone() const = 0;
    virtual reference dereference() const = 0;

  };

  /** @}*/

  // ---- Interface ---------------------

  class PrimClex;
  template<bool>
  class ScelEnumT;
  typedef ScelEnumT<false> ScelEnum;
  class ScelEnumProps;
  class Lattice;
  template<typename T>
  class SupercellEnumerator;

  namespace Completer {
    class EnumOption;
  }

  /** \addtogroup Enumerator
      @{
  */

  /// \brief Base class for generic use of enumerators that may be accessed through the API
  class EnumInterfaceBase {

  public:

    EnumInterfaceBase() {}

    virtual ~EnumInterfaceBase() {}

    virtual std::string help() const = 0;

    virtual std::string name() const = 0;

    virtual int run(PrimClex &primclex, const jsonParser &kwargs, const Completer::EnumOption &enum_opt) const = 0;

    std::unique_ptr<EnumInterfaceBase> clone() const {
      return std::unique_ptr<EnumInterfaceBase>(this->_clone());
    }

  private:

    virtual EnumInterfaceBase *_clone() const = 0;

  };


  /// \brief Used to hold a list of all enumerators that may be accessed via the API
  typedef notstd::unique_cloneable_map<std::string, EnumInterfaceBase> EnumeratorMap;

  /// \brief Use to construct an EnumeratorMap
  inline EnumeratorMap make_enumerator_map() {

    return EnumeratorMap(
    [](const EnumInterfaceBase & e) -> std::string {
      return e.name();
    },
    [](const EnumInterfaceBase & e) -> notstd::cloneable_ptr<EnumInterfaceBase> {
      return notstd::clone(e);
    });
  }

  /// \brief Standardizes parsing casm enum input options to make ScelEnum JSON input
  jsonParser make_enumerator_scel_enum_input(
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt);

  /// \brief Standardizes parsing casm enum input options to make an ScelEnumProps
  ScelEnumProps make_enumerator_scel_enum_props(
    PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt);

  /// \brief Standardizes parsing casm enum input options to make an SupercellEnumerator<Lattice>
  std::unique_ptr<SupercellEnumerator<Lattice> > make_enumerator_superlat_enum(
    PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt);

  /// \brief Standardizes parsing casm enum input options to make an ScelEnum
  std::unique_ptr<ScelEnum> make_enumerator_scel_enum(
    PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt);

  /// \brief Standardizes parsing casm enum filter expressions
  std::vector<std::string> make_enumerator_filter_expr(
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt);

  /// \brief Standardizes insertion from enumerators that construct unique
  /// primitive canonical configurations
  template<typename ScelIterator, typename ConfigEnumConstructor>
  int insert_unique_canon_configs(
    std::string method,
    PrimClex &primclex,
    ScelIterator begin,
    ScelIterator end,
    ConfigEnumConstructor f,
    std::vector<std::string> filter_expr);

  /// \brief Standardizes insertion from enumerators that construct configurations
  template<typename ScelIterator, typename ConfigEnumConstructor>
  int insert_configs(
    std::string method,
    PrimClex &primclex,
    ScelIterator begin,
    ScelIterator end,
    ConfigEnumConstructor f,
    std::vector<std::string> filter_expr,
    bool primitive_only);

  /// \brief Standardizes insertion from enumerators that construct configurations
  template<typename LatticeIterator, typename ConfigEnumConstructor>
  int insert_configs_via_lattice_enum(
    std::string method,
    PrimClex &primclex,
    LatticeIterator begin,
    LatticeIterator end,
    ConfigEnumConstructor f,
    std::vector<std::string> filter_expr,
    bool primitive_only);

  /// \brief Template class to be specialized for each enumerator that may be accessed via the API
  template<typename Derived>
  class EnumInterface : public EnumInterfaceBase {

  public:

    std::string help() const override {
      return Derived::interface_help;
    }

    std::string name() const override {
      return Derived::enumerator_name;
    }

    int run(PrimClex &primclex, const jsonParser &kwargs, const Completer::EnumOption &enum_opt) const override {
      return Derived::run(primclex, kwargs, enum_opt);
    }

    std::unique_ptr<EnumInterfaceBase> clone() const {
      return std::unique_ptr<EnumInterfaceBase>(this->_clone());
    }

  private:

    EnumInterfaceBase *_clone() const override {
      return new EnumInterface<Derived>(*this);
    }


  };

  /** @}*/

}

#endif
