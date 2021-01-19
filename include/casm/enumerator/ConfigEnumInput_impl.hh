#ifndef CASM_ConfigEnumInput_impl
#define CASM_ConfigEnumInput_impl

#include "casm/enumerator/ConfigEnumInput.hh"

namespace CASM {

/** \addtogroup Enumerator
    @{
*/

/// Select by "site_index" or "site_uccoord"
template <typename SiteContainer>
void ConfigEnumInput::select_sites(SiteContainer const &_container) {
  for (auto const &site_index : _container) select_site(site_index);
}

/// Select all sites on multiples sublattices
template <typename SublatticeIndexContainer>
void ConfigEnumInput::select_sublattices(
    SublatticeIndexContainer const &_container) {
  for (auto const &sublattice_index : _container)
    select_sublattice(sublattice_index);
}

/// Returns the subgroup of [group_begin, group_end] that does not cause any
/// permutation between the set of selected and unselected sites of
/// "config_enum_input"
template <typename PermuteIteratorIt>
std::vector<PermuteIterator> make_invariant_subgroup(
    ConfigEnumInput const &config_enum_input, PermuteIteratorIt group_begin,
    PermuteIteratorIt group_end) {
  std::vector<PermuteIterator> invariant_subgroup;
  std::set<Index> const &selected_sites = config_enum_input.sites();

  for (auto it = group_begin; it != group_end; ++it) {
    if (site_indices_are_invariant(*it, selected_sites)) {
      invariant_subgroup.push_back(*it);
    }
  }
  return invariant_subgroup;
}

/** @}*/
}  // namespace CASM

#endif
