#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/clex/Supercell.hh"
#include "casm/enumerator/ConfigEnumInput.hh"

namespace CASM {

  namespace local_impl {
    std::vector<int> max_selected_occupation(ConfigEnumInput const &config_enum_input) {
      auto const &supercell = config_enum_input.configuration().supercell();
      std::vector<int> max_allowed = supercell.max_allowed_occupation();

      std::vector<int> max_allowed_on_selected_sites;
      for(Index i : config_enum_input.sites()) {
        max_allowed_on_selected_sites.push_back(max_allowed[i]);
      }

      return max_allowed_on_selected_sites;
    }

    void set_occupation(Configuration &configuration,
                        std::set<Index> const &site_indices,
                        std::vector<int> const &counter) {
      Index i = 0;
      for(Index site_index : site_indices) {
        configuration.set_occ(site_index, counter[i++]);
      }
    }
  }

  /// Conditionally true for ConfigEnumAllOccupations (true when enumerating on all sites)
  template<>
  bool is_guaranteed_for_database_insert(ConfigEnumAllOccupations const &enumerator) {
    return enumerator.canonical_guarantee();
  }

  /// \brief Construct with a Supercell, using all permutations
  ConfigEnumAllOccupations::ConfigEnumAllOccupations(const ConfigEnumInput &config_enum_input) :
    m_site_index_selection(config_enum_input.sites()),
    m_counter(
      std::vector<int>(config_enum_input.sites().size(), 0),
      local_impl::max_selected_occupation(config_enum_input),
      std::vector<int>(config_enum_input.sites().size(), 1)),
    m_current(notstd::make_cloneable<Configuration>(config_enum_input.configuration())),
    m_enumerate_on_a_subset_of_supercell_sites(
      m_site_index_selection.size() != config_enum_input.configuration().size()) {

    local_impl::set_occupation(*m_current, m_site_index_selection, m_counter);
    reset_properties(*m_current);
    this->_initialize(&(*m_current));

    // Make sure that current() is a primitive canonical config
    if(!_current_is_valid_for_output()) {
      increment();
    }

    // set step to 0
    if(valid()) {
      _set_step(0);
    }
    m_current->set_source(this->source(step()));
  }

  std::string ConfigEnumAllOccupations::name() const {
    return enumerator_name;
  }

  /// Returns true if enumerator is guaranteed to output canonical configurations
  ///
  /// Output is sometimes canonical (when enumerating on all supercell sites) and sometimes not
  /// (when enumerating on a subset of supercell sites). If this returns true, then the output is
  /// guaranteed to be in canonical form.
  bool ConfigEnumAllOccupations::canonical_guarantee() const {
    return !m_enumerate_on_a_subset_of_supercell_sites;
  }

  const std::string ConfigEnumAllOccupations::enumerator_name = "ConfigEnumAllOccupations";

  /// Implements _increment over all occupations
  void ConfigEnumAllOccupations::increment() {

    bool is_valid_config {false};

    while(!is_valid_config && ++m_counter) {
      local_impl::set_occupation(*m_current, m_site_index_selection, m_counter);
      is_valid_config = _current_is_valid_for_output();
    }

    // if while loop exited with valid configuration, increment step number
    if(m_counter.valid()) {
      this->_increment_step();
    }
    // if while loop exited because counter hit end, invalidate this enumerator
    else {
      this->_invalidate();
    }
    m_current->set_source(this->source(step()));
  }

  /// Returns true if current() is primitive and canonical
  bool ConfigEnumAllOccupations::_current_is_valid_for_output() const {
    if(m_enumerate_on_a_subset_of_supercell_sites) {
      return current().is_primitive();
    }
    else {
      return current().is_primitive() && current().is_canonical();
    }
  }

}
