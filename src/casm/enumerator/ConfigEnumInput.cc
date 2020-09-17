#include "casm/clex/Supercell.hh"
#include "casm/clusterography/SupercellClusterOrbits.hh"
#include "casm/enumerator/ConfigEnumInput_impl.hh"
#include "casm/symmetry/PermuteIterator.hh"

namespace CASM {

  /** \addtogroup Enumerator
      @{
  */

  /// Construct with a Configuration and all sites selected
  ConfigEnumInput::ConfigEnumInput(Configuration const &_configuration) :
    m_configuration(_configuration) {
    select_all_sites();
  }

  /// Construct using a Supercell's default Configuration and all sites selected
  ///
  /// Note:
  /// - The Supercell default Configuration is given by `Configuration::zeros(_supercell)`
  ConfigEnumInput::ConfigEnumInput(Supercell const &_supercell) :
    ConfigEnumInput(Configuration::zeros(_supercell)) {}

  /// Construct with a Configuration and specified sites selected
  ///
  /// Note:
  /// - Empty `_site_index_selection` results in no sites being selected
  ConfigEnumInput::ConfigEnumInput(Configuration const &_configuration,
                                   std::set<Index> const &_site_index_selection = {}) :
    m_site_index_selection(_site_index_selection),
    m_configuration(_configuration) {}

  /// Construct using a Supercell's default Configuration and specified sites selected
  ///
  /// Note:
  /// - The Supercell default Configuration is given by `Configuration::zeros(_supercell)`
  /// - Empty `_site_index_selection` results in no sites being selected
  ConfigEnumInput::ConfigEnumInput(Supercell const &_supercell,
                                   std::set<Index> const &_site_index_selection = {}) :
    ConfigEnumInput(Configuration::zeros(_supercell), _site_index_selection) {}

  Configuration const &ConfigEnumInput::configuration() const {
    return m_configuration;
  }

  // /// (deprecated) Prefer using configuration()
  // Configuration const &ConfigEnumInput::config() const {
  //   return configuration();
  // }
  //
  // /// (deprecated) Prefer using configuration()
  // Supercell const &ConfigEnumInput::supercell() const {
  //   return configuration().supercell();
  // }
  //
  // /// (deprecated) Prefer using configuration()
  // ConfigDoF const &ConfigEnumInput::configdof() const {
  //   return configuration().configdof();
  // }

  std::set<Index> const &ConfigEnumInput::sites() const {
    return m_site_index_selection;
  }

  void ConfigEnumInput::clear_sites() {
    m_site_index_selection.clear();
  }

  void ConfigEnumInput::select_all_sites() {
    for(Index i = 0; i < configuration().size(); ++i)
      m_site_index_selection.insert(i);
  }

  void ConfigEnumInput::select_site(Index site_index) {
    m_site_index_selection.insert(site_index);
  }

  void ConfigEnumInput::select_site(UnitCellCoord const &site_uccord) {
    auto const &converter = configuration().supercell().sym_info().unitcellcoord_index_converter();
    m_site_index_selection.insert(converter(site_uccord));
  }

  void ConfigEnumInput::select_sublattice(Index sublattice_index) {
    Index b = sublattice_index;
    Index V = m_configuration.supercell().volume();
    for(Index i = b * V; i < (b + 1)*V; ++i) {
      m_site_index_selection.insert(i);
    }
  }

  /// Returns the subgroup of the configuration factor group that not cause any permutation between
  /// the set of selected and unselected sites of "config_enum_input"
  std::vector<PermuteIterator> make_invariant_group(ConfigEnumInput const &config_enum_input) {

    std::vector<PermuteIterator> invariant_group;
    auto factor_group = config_enum_input.configuration().factor_group();
    std::set<Index> const &selected_sites = config_enum_input.sites();

    for(PermuteIterator const &perm_it : factor_group) {
      if(cluster_site_indices_are_invariant(perm_it, selected_sites)) {
        invariant_group.push_back(perm_it);
      }
    }
    return invariant_group;
  }

  /** @}*/
}
