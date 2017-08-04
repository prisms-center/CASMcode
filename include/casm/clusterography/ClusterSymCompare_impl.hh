#ifndef CASM_ClusterSymCompare_impl
#define CASM_ClusterSymCompare_impl

#include "casm/clusterography/ClusterSymCompare.hh"
#include "casm/crystallography/PrimGrid.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"

namespace CASM {

  /// \brief Return tolerance
  template<typename Derived>
  double ClusterSymCompare<Derived>::tol() const {
    return m_tol;
  }

  /// \brief Constructor
  ///
  /// \param tol Tolerance for invariants_compare of site-to-site distances
  ///
  template<typename Derived>
  ClusterSymCompare<Derived>::ClusterSymCompare(double tol):
    SymCompare<ClusterSymCompare<Derived> >(),
    m_tol(tol) {

  }

  /// \brief Orders 'prepared' elements in the same orbit
  ///
  /// - Returns 'true' to indicate A < B
  /// - Equivalence is indicated by \code !compare(A,B) && !compare(B,A) \endcode
  /// - Assumes elements are 'prepared' before being compared
  /// Implementation:
  /// - First compares by number of sites in cluster
  /// - Then compare all displacements, from longest to shortest
  template<typename Derived>
  bool ClusterSymCompare<Derived>::invariants_compare_impl(const InvariantsType &A, const InvariantsType &B) const {
    return CASM::compare(A, B, tol());
  }

  /// \brief Compares 'prepared' elements
  ///
  /// - Returns 'true' to indicate A < B
  /// - Equivalence is indicated by \code !compare(A,B) && !compare(B,A) \endcode
  /// - Assumes elements are 'prepared' before being compared
  template<typename Derived>
  bool ClusterSymCompare<Derived>::compare_impl(const Element &A, const Element &B) const {
    return A < B;
  }

  /// \brief type-specific way to get position of element
  ///
  /// - Returns traits<Element>::position(el)
  template<typename Derived>
  UnitCellCoord ClusterSymCompare<Derived>::position(const Element &el) {
    return traits<Element>::position(el);
  }


  // -- LocalSymCompare<IntegralCluster> -------------------------------------

  /// \brief Constructor
  ///
  /// \param tol Tolerance for invariants_compare of site-to-site distances
  ///
  template<typename Element>
  AperiodicSymCompare<Element/*, enable_if_integral_position<Element>*/>::
  AperiodicSymCompare(double tol):
    ClusterSymCompare<LocalSymCompare<Element> >(tol) {}

  /// \brief Prepare an element for comparison
  ///
  /// - Returns sorted
  template<typename Element>
  Element AperiodicSymCompare<Element/*, enable_if_integral_position<Element>*/>::
  prepare_impl(Element obj) const {
    return obj.sort();
  }


  // -- PrimPeriodicSymCompare<IntegralCluster> -------------------------------------

  /// \brief Constructor
  ///
  /// \param tol Tolerance for invariants_compare of site-to-site distances
  ///
  template<typename Element>
  PrimPeriodicSymCompare<Element/*, enable_if_integral_position<Element>*/>::
  PrimPeriodicSymCompare(double tol):
    ClusterSymCompare<PrimPeriodicSymCompare<Element> >(tol) {}

  /// \brief Constructor
  ///
  /// \param tol Tolerance for invariants_compare of site-to-site distances
  ///
  template<typename Element>
  PrimPeriodicSymCompare<Element/*, enable_if_integral_position<Element>*/>::
  PrimPeriodicSymCompare(const PrimClex &primclex):
    PrimPeriodicSymCompare(primclex.crystallography_tol()) {}

  /// \brief Prepare an element for comparison
  ///
  /// - Sorts and translates so that obj[0] is in the origin unit cell
  template<typename Element>
  Element PrimPeriodicSymCompare<Element/*, enable_if_integral_position<Element>*/>::
  prepare_impl(Element obj) const {
    if(!obj.size()) {
      return obj;
    }
    obj.sort();
    m_integral_tau = -position(obj).unitcell();
    return obj + m_integral_tau;;
  }


  // -- ScelPeriodicSymCompare<IntegralCluster> -------------------------------------

  /// \brief Constructor
  ///
  /// \param tol Tolerance for invariants_compare of site-to-site distances
  ///
  template<typename Element>
  ScelPeriodicSymCompare<Element/*, enable_if_integral_position<Element>*/>::
  ScelPeriodicSymCompare(const PrimGrid &prim_grid, double tol):
    ClusterSymCompare<ScelPeriodicSymCompare<Element> >(tol),
    m_prim_grid(&prim_grid) {}

  /// \brief Constructor
  ///
  /// \param tol Tolerance for invariants_compare of site-to-site distances
  ///
  template<typename Element>
  ScelPeriodicSymCompare<Element/*, enable_if_integral_position<Element>*/>::
  ScelPeriodicSymCompare(const Supercell &scel):
    ScelPeriodicSymCompare(scel.prim_grid(), scel.crystallography_tol()) {}

  /// \brief Prepare an element for comparison
  ///
  /// - Sorts UnitCellCoord and translates so that obj[0] is within the supercell
  template<typename Element>
  Element ScelPeriodicSymCompare<Element/*, enable_if_integral_position<Element>*/>::
  prepare_impl(Element obj) const {
    if(!obj.size()) {
      return obj;
    }
    obj.sort();
    auto pos = position(obj);
    m_integral_tau = pos.unitcell() - m_prim_grid->within(pos).unitcell();
    return obj + m_integral_tau;
  }

}

#endif
