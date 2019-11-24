#ifndef CASM_SubOrbits_impl
#define CASM_SubOrbits_impl

#include "casm/symmetry/InvariantSubgroup_impl.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SubOrbits.hh"

namespace CASM {

  // --- MakeSubOrbitGenerators ---
  //

  template<typename CopyApplyType>
  template<typename OrbitType, typename ElementOutputIterator>
  ElementOutputIterator MakeSubOrbitGenerators<CopyApplyType>::operator()(
    const OrbitType &orbit,
    ElementOutputIterator result) const {

    SymGroup invariant_subgroup = make_invariant_subgroup(orbit);
    return (*this)(orbit.prototype(), invariant_subgroup, result);
  }

  template<typename CopyApplyType>
  template<typename Element, typename SymCompareType, typename ElementOutputIterator>
  ElementOutputIterator MakeSubOrbitGenerators<CopyApplyType>::operator()(
    const Element &element,
    const SymCompareType &sym_compare,
    ElementOutputIterator result) const {

    SymGroup invariant_subgroup = make_invariant_subgroup(element, m_group, sym_compare);
    return (*this)(element, invariant_subgroup, result);
  }

  template<typename CopyApplyType>
  template<typename Element, typename ElementOutputIterator>
  ElementOutputIterator MakeSubOrbitGenerators<CopyApplyType>::operator()(
    const Element &element,
    const SymGroup &invariant_subgroup,
    ElementOutputIterator result) const {

    // find "max" SymOp in each coset of the subgroup,
    //   excluding those cosets that cause duplicate generating elements
    for(Index i = 0; i < m_group.size(); ++i) {
      const SymOp &test_op = m_group[i];
      auto lambda = [&](const SymOp & op) {
        for(const auto &el_op : invariant_subgroup) {
          if(test_op.index() < (op * test_op * el_op).index()) {
            return true;
          }
        }
        return false;
      };
      // if test_op is max
      if(std::none_of(m_subgroup.begin(), m_subgroup.end(), lambda)) {
        // apply to element and construct suborbit generator
        *result++ = this->m_copy_apply_f(test_op, element);
      }
    }
    return result;
  }

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to group -> subgroup symmetry breaking
  template<typename Element, typename CopyApplyElementType, typename ElementOutputIterator>
  ElementOutputIterator make_suborbit_generators(
    const Element &element,
    const SymGroup &invariant_subgroup,
    const SymGroup &group,
    const SymGroup &subgroup,
    const CopyApplyElementType &copy_apply_f,
    ElementOutputIterator result) {

    MakeSubOrbitGenerators<CopyApplyElementType> f{group, subgroup, copy_apply_f};
    return f(element, invariant_subgroup, result);
  }

  /// \brief Output the orbit generators necessary to construct the sub-orbits
  /// corresponding to group -> subgroup symmetry breaking
  template<typename OrbitType, typename CopyApplyElementType, typename ElementOutputIterator>
  ElementOutputIterator make_suborbit_generators(
    const OrbitType &orbit,
    const SymGroup &group,
    const SymGroup &subgroup,
    const CopyApplyElementType &copy_apply_f,
    ElementOutputIterator result) {

    MakeSubOrbitGenerators<CopyApplyElementType> f{group, subgroup, copy_apply_f};
    return f(orbit, result);
  }

}

#endif
