#ifndef CASM_InvariantSubgroup_impl
#define CASM_InvariantSubgroup_impl

#include "casm/symmetry/InvariantSubgroup.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymCompare.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/clusterography/ClusterSymCompare.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {

  using xtal::Coordinate;

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
  ///   if(sym_compare.equal(e, sym_compare.prepare(sym_compare.copy_apply(op, e)))) {
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
      if(sym_compare.equal(e, sym_compare.prepare(sym_compare.copy_apply(op, e)))) {
        result.push_back(sym_compare.spatial_transform()*op);
      }
    }
    if(result[0].index() != 0) {
      throw std::runtime_error("Error in make_invariant_subgroup (0): First element is not identity.");
    }
    return result;
  }

  /// \brief Construct the subgroup that leaves an element of the orbit unchanged
  ///
  /// Does not include translations as determined from 'sym_compare'
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
    const auto &map = orbit.equivalence_map();
    result.set_lattice(
      map[0][0].master_group().lattice());
    for(Index i = 0; i < orbit.equivalence_map()[0].size(); ++i) {
      result.push_back(map[element_index][0]*map[0][i]*map[element_index][0].inverse());
    }
    if(result[0].index() != 0) {
      throw std::runtime_error("Error in make_invariant_subgroup (1): First element is not identity.");
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
             scel.sym_info().permute_begin(),
             scel.sym_info().permute_end());
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
      scel.primclex().shared_prim(),
      scel.transf_mat(),
      scel.crystallography_tol());
    Element e(sym_compare.prepare(element));
    std::vector<PermuteIterator> result;
    Coordinate coord(scel.prim().lattice());
    auto it = begin;
    while(it != end) {
      auto test = sym_compare.prepare(sym_compare.copy_apply(it.sym_op(), e));
      if(sym_compare.equal(test, e)) {
        coord.cart() = sym_compare.spatial_transform().integral_tau();
        auto trans_it = scel.sym_info().permute_it(0, xtal::UnitCell::from_coordinate(coord));
        result.push_back(trans_it * it);
      }
      ++it;
    }
    if(result[0].sym_op().index() != 0) {
      throw std::runtime_error("Error in make_invariant_subgroup (2): First element is not identity.");
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
      scel.primclex().shared_prim(),
      scel.transf_mat(),
      scel.crystallography_tol());
    Element e(sym_compare.prepare(element));
    std::vector<PermuteIterator> result;
    Coordinate coord(scel.prim().lattice());
    auto it = begin;
    while(it != end) {
      auto test = sym_compare.prepare(sym_compare.copy_apply(it->sym_op(), e));
      if(sym_compare.equal(test, e)) {
        coord.cart() = sym_compare.spatial_transform().integral_tau();
        auto trans_it = scel.sym_info().permute_it(0, xtal::UnitCell::from_coordinate(coord));
        result.push_back(trans_it * (*it));
      }
      ++it;
    }
    if(result[0].sym_op().index() != 0) {
      throw std::runtime_error("Error in make_invariant_subgroup (3): First element is not identity.");
    }
    return result;
  }

  /// \brief Construct the subgroup of permutations that leaves a Supercell unchanged
  ///
  /// \param scel_A Supercell find subgroup of permutations that leave scel_A unchanged
  /// \param scel_B Supercell associated with the supergroup [begin, end)
  /// \param begin,end Range of PermuteIterator describing the supergroup
  template<typename PermuteIteratorIt>
  std::vector<PermuteIterator> make_invariant_subgroup(
    const Supercell &scel_A,
    const Supercell &scel_B,
    PermuteIteratorIt begin,
    PermuteIteratorIt end) {

    const SymGroup &scel_A_fg = scel_A.factor_group();
    const SymGroup &scel_B_fg = scel_B.factor_group();

    auto find_fg_op = [&](const PermuteIterator & scel_B_it) {
      Index master_fg_index = scel_B_fg[scel_B_it.factor_group_index()].index();
      return std::any_of(
               scel_A_fg.begin(),
               scel_A_fg.end(),
      [&](const SymOp & op) {
        return op.index() == master_fg_index;
      });
    };

    std::vector<PermuteIterator> subgroup;
    std::copy_if(begin, end, std::back_inserter(subgroup), find_fg_op);
    if(subgroup[0].sym_op().index() != 0) {
      throw std::runtime_error("Error in make_invariant_subgroup (4): First element is not identity.");
    }
    return subgroup;
  }

}

#endif
