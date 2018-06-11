#ifndef CASM_ClusterSymCompare_impl
#define CASM_ClusterSymCompare_impl

#include "casm/clusterography/ClusterSymCompare.hh"
#include "casm/symmetry/SymPermutation.hh"
#include "casm/crystallography/PrimGrid.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"

namespace CASM {

  /// \brief Return tolerance
  template<typename Base>
  double ClusterSymCompare<Base>::tol() const {
    return m_tol;
  }

  /// \brief Constructor
  ///
  /// \param tol Tolerance for invariants_compare of site-to-site distances
  ///
  template<typename Base>
  ClusterSymCompare<Base>::ClusterSymCompare(double tol):
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
  template<typename Base>
  bool ClusterSymCompare<Base>::invariants_compare_impl(const InvariantsType &A, const InvariantsType &B) const {
    return CASM::compare(A, B, tol());
  }

  /// \brief Compares 'prepared' elements
  ///
  /// - Returns 'true' to indicate A < B
  /// - Equivalence is indicated by \code !compare(A,B) && !compare(B,A) \endcode
  /// - Assumes elements are 'prepared' before being compared
  template<typename Base>
  bool ClusterSymCompare<Base>::compare_impl(const Element &A, const Element &B) const {
    return A < B;
  }


  /// \brief Returns transformation that takes 'obj' to its prepared (canonical) form
  ///
  /// - For now returns pointer to SymPermutation object that encodes permutation due to sorting elements
  template<typename Base>
  std::unique_ptr<SymOpRepresentation> ClusterSymCompare<Base>::canonical_transform_impl(Element const &obj)const {
    return std::unique_ptr<SymOpRepresentation>(new SymPermutation(obj.sort_permutation()));
  }

  /// \brief type-specific way to get position of element
  ///
  /// - Returns traits<Element>::position(el)
  template<typename Base>
  UnitCellCoord ClusterSymCompare<Base>::position(const Element &el) {
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
    ClusterSymCompare<SymCompare<CRTPBase<AperiodicSymCompare<Element>>>>(tol) {}

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
    ClusterSymCompare<SymCompare<CRTPBase<PrimPeriodicSymCompare<Element>>>>(tol) {}

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
    this->m_integral_tau = -position(obj).unitcell();
    return obj + this->m_integral_tau;;
  }


  // -- ScelPeriodicSymCompare<IntegralCluster> -------------------------------------

  /// \brief Constructor
  ///
  /// \param tol Tolerance for invariants_compare of site-to-site distances
  ///
  template<typename Element>
  ScelPeriodicSymCompare<Element/*, enable_if_integral_position<Element>*/>::
  ScelPeriodicSymCompare(const PrimGrid &prim_grid, double tol):
    ClusterSymCompare<SymCompare<CRTPBase<ScelPeriodicSymCompare<Element>>>>(tol),
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
    this->m_integral_tau = m_prim_grid->within(pos).unitcell() - pos.unitcell();
    return obj + this->m_integral_tau;
  }


  // -- WithinScelSymCompare<IntegralCluster> -------------------------------------

  /// \brief Constructor
  ///
  /// \param tol Tolerance for invariants_compare of site-to-site distances
  ///
  template<typename Element>
  WithinScelSymCompare<Element/*, enable_if_integral_position<Element>*/>::
  WithinScelSymCompare(const PrimGrid &prim_grid, double tol):
    ClusterSymCompare<SymCompare<CRTPBase<WithinScelSymCompare<Element>>>>(tol),
    m_prim_grid(&prim_grid) {}

  /// \brief Constructor
  ///
  /// \param tol Tolerance for invariants_compare of site-to-site distances
  ///
  template<typename Element>
  WithinScelSymCompare<Element/*, enable_if_integral_position<Element>*/>::
  WithinScelSymCompare(const Supercell &scel):
    WithinScelSymCompare<Element>(scel.prim_grid(), scel.crystallography_tol()) {}

  /// \brief Prepare an element for comparison
  ///
  /// - Puts all sites within the supercell, then sorts
  template<typename Element>
  Element WithinScelSymCompare<Element/*, enable_if_integral_position<Element>*/>::
  prepare_impl(Element obj) const {
    if(!obj.size()) {
      return obj;
    }
    for(Index i = 0; i < obj.size(); ++i) {
      obj[i] = m_prim_grid->within(obj[i]);
    }
    obj.sort();
    return obj;
  }
}

#endif
