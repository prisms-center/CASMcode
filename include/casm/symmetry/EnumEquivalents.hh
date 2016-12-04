#ifndef CASM_EnumEquivalents
#define CASM_EnumEquivalents

#include "casm/container/InputEnumerator.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

  /** \defgroup EnumEquivalents
   *
   *  \brief Enumerate equivalent objects
   *
   *  \ingroup Symmetry
   *  \ingroup Enumerator
   *
   *  @{
  */

  template <
    typename ObjectType,
    typename SymOpIterator,
    typename SymOpType,
    typename SymOpCompare >
  class EnumEquivalents;

  namespace CASM_TMP {

    template <
      typename ObjectType,
      typename SymOpIterator,
      typename SymOpType,
      typename SymOpCompare >
    struct traits<EnumEquivalents<ObjectType, SymOpIterator, SymOpType, SymOpCompare> > {
      static const std::string name;
    };

  }

  /// \brief Enumerate over all equivalents
  ///
  /// Assumes the following exist:
  /// - ObjectType copy_apply(const SymOpType&, const ObjectType&)
  /// - bool operator<(SymOpType& A, SymOpType& B);
  ///
  ///
  /// \ingroup EnumEquivalents
  ///
  template <
    typename ObjectType,
    typename SymOpIterator,
    typename SymOpType = typename std::iterator_traits<SymOpIterator>::value_type,
    typename SymOpCompare = std::less<SymOpType> >
  class EnumEquivalents : public InputEnumeratorBase<ObjectType> {

    // -- Required members -------------------

  public:

    /// \brief Constructor
    ///
    /// \param generating_obj The object to enumerate equivalents of
    /// \param begin,end A range of SymOpType
    /// \param f A functor that outputs the subgroup of [begin,end) that leaves
    ///          the generating_obj invariant
    ///
    /// Expected form for f:
    /// \code
    /// template<typename SymOpIterator, typename SymOpOutputIterator>
    /// SymOpOutputIterator f(
    ///   const ObjectType &generating_obj,
    ///   SymOpIterator begin,
    ///   SymOpIterator end,
    ///   SymOpOutputIterator result);
    /// \endcode
    ///
    template<typename MakeInvariantSubgroup>
    EnumEquivalents(const ObjectType &generating_obj,
                    SymOpIterator begin,
                    SymOpIterator end,
                    MakeInvariantSubgroup f = MakeInvariantSubgroup(),
                    SymOpCompare compare = SymOpCompare()) :
      m_initial(notstd::clone(generating_obj)),
      m_current(m_initial),
      m_symop_it(begin),
      m_symop_end(end),
      m_compare(compare) {

      f(generating_obj, begin, end, std::back_inserter(m_invariant_subgroup));

      this->_initialize(&(*m_current));

      if(!_check(*m_symop_it)) {
        increment();
      }

      // set step to 0
      if(this->valid()) {
        this->_set_step(0);
      }
    }


  public:

    // -- Unique -------------------

    const ObjectType &generator() const {
      return *m_initial;
    }

    const SymOpType &sym_op() const {
      return *m_symop_it;
    }

    const std::vector<SymOpType> &invariant_subgroup() const {
      return m_invariant_subgroup;
    }

  private:

    // -- Required members -------------------

    /// Implement increment
    void increment() override {

      while(++m_symop_it != m_symop_end && !_check(*m_symop_it)) {
        // increment
      }
      if(m_symop_it == m_symop_end) {
        this->_invalidate();
      }
      else {

        // generate equivalent by applying symop
        *m_current = copy_apply(*m_symop_it, *m_initial);

        // increment
        this->_increment_step();
      }
    }

    // -- Unique -------------------

    /// Returns true for one element (the "minimum" SymOp) of the left
    /// cosets of the invariant subgroup
    bool _check(const SymOpType &X) const {
      return std::none_of(
               invariant_subgroup().begin(),
               invariant_subgroup().end(),
      [&](const SymOpType & Bi) {
        return m_compare(X * Bi, X);
      });
    }

    notstd::cloneable_ptr<ObjectType> m_initial;
    notstd::cloneable_ptr<ObjectType> m_current;
    SymOpIterator m_symop_it;
    SymOpIterator m_symop_end;
    std::vector<SymOpType> m_invariant_subgroup;
    SymOpCompare m_compare;
  };

  namespace CASM_TMP {

    template <
      typename ObjectType,
      typename SymOpIterator,
      typename SymOpType,
      typename SymOpCompare >
    const std::string traits<EnumEquivalents<ObjectType, SymOpIterator, SymOpType, SymOpCompare> >::name = "EnumEquivalents";

  }

}

#endif
