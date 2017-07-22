#ifndef CASM_InvariantSubgroup_impl
#define CASM_InvariantSubgroup_impl

#include "casm/symmetry/InvariantSubgroup.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymCompare.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/clex/Supercell.hh"

namespace CASM {

  /// \brief Construct the subgroup that leaves an element unchanged
  ///
  /// Includes translations as determined from 'sym_compare'
  ///
  /// Implementation:
  /// \code
  /// Element e(sym_compare.prepare(element));
  /// SymGroup result = generating_grp;
  /// result.clear();
  /// for(const auto &op : generating_grp) {
  ///   if(sym_compare.equal(e, sym_compare.prepare(copy_apply(op, e)))) {
  ///    result.push_back(sym_compare.translation(element.prim())*op);
  ///   }
  /// }
  /// return result;
  /// \endcode
  template<typename Element, typename SymCompareType>
  SymGroup make_invariant_subgroup(const Element &element,
                                   const SymGroup &generating_grp,
                                   const SymCompareType &sym_compare) {
    Element e(sym_compare.prepare(element));
    SymGroup result = generating_grp;
    result.clear();
    for(const auto &op : generating_grp) {
      if(sym_compare.equal(e, sym_compare.prepare(copy_apply(op, e)))) {
        result.push_back(sym_compare.translation(element.prim())*op);
      }
    }
    return result;
  }

  /// \brief Construct the subgroup that leaves an element of the orbit unchanged
  ///
  /// The equivalence map is:
  ///   element(i) compares equivalent to prototype().copy_apply(equivalence_map[i][j]) for all j
  ///
  /// The subgroup that leaves element 'i' of the orbit unchanged is, for all j:
  ///   eq_map[i][0]*eq_map[0][j]*eq_map[i][0].inverse()
  ///
  template<typename OrbitType>
  SymGroup make_invariant_subgroup(const OrbitType &orbit, Index element_index) {
    SymGroup result;
    result.set_lattice(orbit.prototype().lattice());
    const auto &map = orbit.equivalence_map();
    for(Index i = 0; i < orbit.equivalence_map()[0].size(); ++i) {
      result.push_back(map[element_index][0]*map[0][i]*map[element_index][0].inverse());
    }
    return result;
  }

  /// \brief Construct the subgroup of permutations that leaves an element unchanged
  ///
  /// Uses comparison defined by (and include translation):
  /// \code
  /// ScelPeriodicSymCompare<Element> sym_compare(
  ///    scel.prim_grid(),
  ///    scel.crystallography_tol());
  /// \endcode
  ///
  template<typename Element>
  std::vector<PermuteIterator> make_invariant_subgroup(const Element &element, const Supercell &scel) {
    return make_invariant_subgroup(
             element,
             scel,
             scel.permute_begin(),
             scel.permute_end());
  }

  /// \brief Construct the subgroup of permutations that leaves an element unchanged
  ///
  /// Uses comparison defined by (and include translation):
  /// \code
  /// ScelPeriodicSymCompare<Element> sym_compare(
  ///    scel.prim_grid(),
  ///    scel.crystallography_tol());
  /// \endcode
  ///
  template<typename Element>
  std::vector<PermuteIterator> make_invariant_subgroup(
    const Element &element,
    const Supercell &scel,
    PermuteIterator begin,
    PermuteIterator end) {

    ScelPeriodicSymCompare<Element> sym_compare(
      scel.prim_grid(),
      scel.crystallography_tol());
    Element e(sym_compare.prepare(element));
    std::vector<PermuteIterator> result;
    auto it = begin;
    while(it != end) {
      auto test = sym_compare.prepare(copy_apply(it.sym_op(), e));
      if(sym_compare.equal(test, e)) {
        auto trans_it = scel.permute_it(0, scel.prim_grid().find(sym_compare.integral_tau()));
        result.push_back(trans_it * it);
      }
      ++it;
    }
    return result;
  }

  /// \brief Construct the subgroup of permutations that leaves an element unchanged
  ///
  /// Uses comparison defined by (and include translation):
  /// \code
  /// ScelPeriodicSymCompare<Element> sym_compare(
  ///    scel.prim_grid(),
  ///    scel.crystallography_tol());
  /// \endcode
  ///
  template<typename Element, typename PermuteIteratorIt>
  std::vector<PermuteIterator> make_invariant_subgroup(
    const Element &element,
    const Supercell &scel,
    PermuteIteratorIt begin,
    PermuteIteratorIt end) {

    ScelPeriodicSymCompare<Element> sym_compare(
      scel.prim_grid(),
      scel.crystallography_tol());
    Element e(sym_compare.prepare(element));
    std::vector<PermuteIterator> result;
    auto it = begin;
    while(it != end) {
      auto test = sym_compare.prepare(copy_apply(it->sym_op(), e));
      if(sym_compare.equal(test, e)) {
        auto trans_it = scel.permute_it(0, scel.prim_grid().find(sym_compare.integral_tau()));
        result.push_back(trans_it * (*it));
      }
      ++it;
    }
    return result;
  }

}

#endif
