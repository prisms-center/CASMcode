#include "casm/clex/ConfigEnumAllOccupations.hh"

#include "casm/clex/Supercell.hh"
#include "casm/enumerator/ConfigEnumInput.hh"

namespace CASM {

namespace local_impl {
std::vector<int> max_selected_occupation(
    ConfigEnumInput const &config_enum_input) {
  auto const &supercell = config_enum_input.configuration().supercell();
  Eigen::VectorXi max_allowed = supercell.max_allowed_occupation();

  std::vector<int> max_allowed_on_selected_sites;
  for (Index i : config_enum_input.sites()) {
    max_allowed_on_selected_sites.push_back(max_allowed[i]);
  }

  return max_allowed_on_selected_sites;
}

void set_occupation(Configuration &configuration,
                    std::set<Index> const &site_indices,
                    std::vector<int> const &counter) {
  Index i = 0;
  for (Index site_index : site_indices) {
    configuration.set_occ(site_index, counter[i++]);
  }
}
}  // namespace local_impl

/// \brief Conditionally true for ConfigEnumAllOccupations
template <>
bool is_guaranteed_for_database_insert(
    ConfigEnumAllOccupations const &enumerator) {
  return enumerator.primitive_canonical_guarantee();
}

/// \brief Construct with a ConfigEnumInput, specifying which sites
///     to enumerate on and which to keep fixed
///
/// Note:
/// - The output configurations are set to fixed occupation matching
///   `config_enum_input.configuration()` on all sites not included in
///   `config_enum_input.sites()`.
/// - All allowed occupation values are enumerated on the selected sites
///   (`config_enum_input.sites()`), but Configuration are only output if the
///   are primitive.
/// - If all sites are selected, then only primitive and canonical
///   Configuration are output. This can be checked with
///   `this->canonical_guarantee()`.
/// - Otherwise, if a subset of sites is selected only primitive
///   Configuration are output
/// - Use the alternative constructor for direct control over whether
///   output configurations are primitive and/or canonical
ConfigEnumAllOccupations::ConfigEnumAllOccupations(
    const ConfigEnumInput &config_enum_input)
    : m_site_index_selection(config_enum_input.sites()),
      m_counter(std::vector<int>(config_enum_input.sites().size(), 0),
                local_impl::max_selected_occupation(config_enum_input),
                std::vector<int>(config_enum_input.sites().size(), 1)),
      m_current(notstd::make_cloneable<Configuration>(
          config_enum_input.configuration())),
      m_enumerate_on_a_subset_of_supercell_sites(
          m_site_index_selection.size() !=
          config_enum_input.configuration().size()) {
  if (m_enumerate_on_a_subset_of_supercell_sites) {
    m_primitive_only = true;
    m_canonical_only = false;
  } else {
    m_primitive_only = true;
    m_canonical_only = true;
  }
  local_impl::set_occupation(*m_current, m_site_index_selection, m_counter);
  reset_properties(*m_current);
  this->_initialize(&(*m_current));

  // Make sure that current() is a primitive canonical config
  if (!_current_is_valid_for_output()) {
    increment();
  }

  // set step to 0
  if (valid()) {
    _set_step(0);
  }
  m_current->set_source(this->source(step()));
}

/// \brief Constructor allowing direct control of whether
///     non-primitive and non-canonical Configuration are enumerated
///
/// \param config_enum_input Specifies the background configuration
///     and sites where all allowed occupation values are enumerated
/// \param primitive_only If true, only primitive configuration are
///     output
/// \param canonical_only If true, only canonical configuration are
///     output
///
/// Note:
/// - The output configurations are set to fixed occupation matching
///   `config_enum_input.configuration()` on all sites not included in
///   `config_enum_input.sites()`.
/// - All allowed occupation values are enumerated on the selected sites
///   (`config_enum_input.sites()`)
/// - If `primitive_only==true`, then non-primitive configurations are
///   skipped.
/// - If `canonical_only==true`, then non-canonical configurations are
///   skipped.
ConfigEnumAllOccupations::ConfigEnumAllOccupations(
    const ConfigEnumInput &config_enum_input, bool primitive_only,
    bool canonical_only)
    : m_site_index_selection(config_enum_input.sites()),
      m_counter(std::vector<int>(config_enum_input.sites().size(), 0),
                local_impl::max_selected_occupation(config_enum_input),
                std::vector<int>(config_enum_input.sites().size(), 1)),
      m_current(notstd::make_cloneable<Configuration>(
          config_enum_input.configuration())),
      m_enumerate_on_a_subset_of_supercell_sites(
          m_site_index_selection.size() !=
          config_enum_input.configuration().size()),
      m_primitive_only(primitive_only),
      m_canonical_only(canonical_only) {
  local_impl::set_occupation(*m_current, m_site_index_selection, m_counter);
  reset_properties(*m_current);
  this->_initialize(&(*m_current));

  // Make sure that current() is a primitive canonical config
  if (!_current_is_valid_for_output()) {
    increment();
  }

  // set step to 0
  if (valid()) {
    _set_step(0);
  }
  m_current->set_source(this->source(step()));
}

std::string ConfigEnumAllOccupations::name() const { return enumerator_name; }

/// \brief Returns true if enumerator is guaranteed to output
///     primitive & canonical configurations only
bool ConfigEnumAllOccupations::primitive_canonical_guarantee() const {
  return m_primitive_only && m_canonical_only;
}

const std::string ConfigEnumAllOccupations::enumerator_name =
    "ConfigEnumAllOccupations";

/// Implements _increment over all occupations
void ConfigEnumAllOccupations::increment() {
  bool is_valid_config{false};

  while (!is_valid_config && ++m_counter) {
    local_impl::set_occupation(*m_current, m_site_index_selection, m_counter);
    is_valid_config = _current_is_valid_for_output();
  }

  // if while loop exited with valid configuration, increment step number
  if (m_counter.valid()) {
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
  if (m_primitive_only && !current().is_primitive()) {
    return false;
  }
  if (m_canonical_only && !current().is_canonical()) {
    return false;
  }
  return true;
}

}  // namespace CASM
