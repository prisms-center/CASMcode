#ifndef CASM_SymCompare
#define CASM_SymCompare

#include <utility>
#include <memory>
#include "casm/CASM_global_definitions.hh"

namespace CASM {

  /* -- SymCompare Declarations --------------------------- */

  class SymOp;

  /// \brief Abstract base class for implementing symmetry properties and
  /// canonical form determination
  template<typename Element>
  class SymCompare {

  public:

    /// \brief Prepare an element for comparison
    ///
    /// - For instance, sort and/or translate a cluster so comparison may be
    ///   performed more efficiently
    virtual Element prepare(Element obj) const = 0;

    /// \brief Orders 'prepared' elements in the same orbit
    ///
    /// - Returns 'true' to indicate A < B
    /// - Assumes elements are 'prepared' before being compared
    virtual bool intra_orbit_compare(const Element &A, const Element &B) const = 0;

    /// \brief Check equivalence of elements in the same orbit
    ///
    /// \returns \code !intra_orbit_compare(A,B) && !intra_orbit_compare(B,A) \endcode
    ///
    /// - Assumes elements are 'prepared' before being compared
    bool intra_orbit_equal(const Element &A, const Element &B) const {
      return !intra_orbit_compare(A, B) && !intra_orbit_compare(B, A);
    }

    /// \brief Orders orbit prototypes in canonical form
    ///
    /// - Returns 'true' to indicate A < B
    /// - Equivalence is indicated by \code !compare(A,B) && !compare(B,A) \endcode
    /// - Assumes elements are in canonical form
    virtual bool inter_orbit_compare(const Element &A, const Element &B) const = 0;

    /// \brief Check equivalence of prototypes in different orbit
    ///
    /// \returns \code !compare(A,B) && !compare(B,A) \endcode
    ///
    /// - Assumes elements are in canonical form
    bool inter_orbit_equal(const Element &A, const Element &B) const {
      return !inter_orbit_compare(A, B) && !inter_orbit_compare(B, A);
    }

    /// \brief Apply symmetry to this
    ///
    /// - For instance, this might change the direction of periodicity, affecting
    ///   the outcome of 'prepare' and/or 'compare'
    SymCompare<Element> &apply_sym(const SymOp &op) {
      this->_apply_sym(op);
      return *this;
    }

    /// \brief Public non-virtual clone
    std::unique_ptr<SymCompare<Element> > clone() const {
      return std::unique_ptr<SymCompare<Element> >(this->_clone());
    }

  private:

    /// \brief Private virtual apply_sym
    virtual void _apply_sym(const SymOp &op) const = 0;

    /// \brief Private virtual clone
    virtual SymCompare<Element> *_clone() const = 0;

  };

  /// \brief Template class to be specialized for comparisons with aperiodic symmetry
  template<typename Element>
  class LocalSymCompare {};

  /// \brief Template class to be specialized for comparisons with periodic symmetry
  /// of the primitive lattice
  template<typename Element>
  class PrimPeriodicSymCompare {};

  /// \brief Template class to be specialized for comparisons with periodic symmetry
  /// of the supercell lattice
  template<typename Element>
  class ScelPeriodicSymCompare {};


  /// \brief Return subgroup that leaves an element unchanged
  ///
  /// All SymOp such that:
  /// \code
  /// Element e = sym_compare.prepare(generating_element);
  /// Element test = sym_compare.prepare(copy_apply(op, e));
  /// sym_compare.equal(e, test) == true
  /// \endcode
  template<typename Element>
  SymGroup invariant_subgroup(const Element &element,
                              const SymGroup &generating_grp,
                              const SymCompare<Element> &sym_compare) {
    Element e(sym_compare.prepare(element));
    SymGroup result;
    for(const auto &op : generating_grp) {
      if(sym_compare.intra_orbit_equal(e, sym_compare.prepare(copy_apply(op, e)))) {
        result.push_back(op);
      }
    }
    return result;
  }

}

#endif
