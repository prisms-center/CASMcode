#ifndef CASM_ConfigEnumInput
#define CASM_ConfigEnumInput

#include <vector>
#include <set>
#include <string>

#include "casm/global/definitions.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {

  /** \addtogroup Enumerator
      @{
  */

  /// Stores a Configuration and a selection of sites in the Configuration
  ///
  /// Site naming convention:
  /// - "sites" and "site_index": Linear index of sites in ConfigDoF
  /// - "site_uccoord": UnitCellCoord index of a site in the Configuration
  /// - "sublattice_index": Linear index of basis site in the prim Structure
  class ConfigEnumInput {
  public:

    /// Construct with a Configuration and all sites selected
    ConfigEnumInput(Configuration const &_configuration);

    /// Construct using a Supercell's default Configuration and all sites selected
    ConfigEnumInput(Supercell const &_supercell);

    /// Construct with a Configuration and initial site selection
    ConfigEnumInput(Configuration const &_configuration,
                    std::set<Index> const &_site_index_selection);

    /// Construct using a Supercell's default Configuration and initial site selection
    ConfigEnumInput(Supercell const &_supercell,
                    std::set<Index> const &_site_index_selection);

    Configuration const &configuration() const;

    // /// (deprecated) Prefer using configuration()
    // Configuration const &config() const;
    //
    // /// (deprecated) Prefer using configuration().supercell()
    // Supercell const &supercell() const;
    //
    // /// (deprecated) Prefer using configuration().configdof()
    // ConfigDoF const &configdof() const;

    /// Const access selected sites, speficied by "site_index" (linear index into ConfigDoF)
    std::set<Index> const &sites() const;

    /// Clear site selection. After condition is no sites are selected.
    void clear_sites();

    void select_all_sites();

    void select_site(Index site_index);

    void select_site(UnitCellCoord const &site_uccord);

    /// Select sites by "site_index" or "site_uccoord"
    template<typename SiteContainer>
    void select_sites(SiteContainer const &_container);

    /// Select all sites on a sublattice
    void select_sublattice(Index sublattice_index);

    /// Select all sites on multiples sublattices
    template<typename SublatticeIndexContainer>
    void select_sublattices(SublatticeIndexContainer const &_container);

  private:

    std::set<Index> m_site_index_selection;
    Configuration m_configuration;
  };

  /// Returns, as a std::vector<PermuteIterator>, the subgroup of
  /// `config_enum_input.configuraiton().factor_group()` that does cause any permutation between the
  /// set of selected and unselected sites of "config_enum_input"
  std::vector<PermuteIterator> make_invariant_group(ConfigEnumInput const &config_enum_input);

  /** @}*/
}

#endif
