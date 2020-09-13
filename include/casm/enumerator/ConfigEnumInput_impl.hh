#ifndef CASM_ConfigEnumInput_impl
#define CASM_ConfigEnumInput_impl

#include "casm/enumerator/ConfigEnumInput.hh"

namespace CASM {

  /** \addtogroup Enumerator
      @{
  */

  /// Select by "site_index" or "site_uccoord"
  template<typename SiteContainer>
  void ConfigEnumInput::select_sites(SiteContainer const &_container) {
    for(auto const &site_index : _container)
      select_site(site_index);
  }

  /// Select all sites on multiples sublattices
  template<typename SublatticeIndexContainer>
  void ConfigEnumInput::select_sublattices(SublatticeIndexContainer const &_container) {
    for(auto const &sublattice_index : _container)
      select_sublattice(sublattice_index);
  }

  /** @}*/
}

#endif
