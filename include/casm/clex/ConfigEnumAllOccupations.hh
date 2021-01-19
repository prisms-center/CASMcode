#ifndef CASM_ConfigEnumAllOccupations
#define CASM_ConfigEnumAllOccupations

#include "casm/clex/Configuration.hh"
#include "casm/container/Counter.hh"
#include "casm/enumerator/InputEnumerator.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

class ConfigEnumAllOccupations;
class ConfigEnumInput;

/** \defgroup ConfigEnumGroup Configuration Enumerators
 *  \ingroup Configuration
 *  \ingroup Enumerator
 *  \brief Enumerates Configuration
 *  @{
 */

/// Conditionally true for ConfigEnumAllOccupations (true when enumerating on
/// all sites)
template <>
bool is_guaranteed_for_database_insert(
    ConfigEnumAllOccupations const &enumerator);

/// Enumerate over all possible occupations on particular sites in a
/// Configuration
///
class ConfigEnumAllOccupations : public InputEnumeratorBase<Configuration> {
  // -- Required members -------------------

 public:
  /// Construct with a ConfigEnumInput, specifying which sites to enumerate on
  /// and which to keep fixed
  ///
  /// Note:
  /// - The output configurations are set to fixed occupation matching
  ///   `config_enum_input.configuration()` on all sites not included in
  ///   `config_enum_input.sites()`.
  /// - All allowed occupation values are enumerated on the selected sites
  ///   (`config_enum_input.sites()`), but Configuration are only output if the
  ///   are primitive.
  /// - If all sites are selected, then only canonical Configuration are output.
  /// This can be
  ///   checked with `this->canonical_guarantee()`.
  ConfigEnumAllOccupations(ConfigEnumInput const &config_enum_input);

  std::string name() const override;

  /// Returns true if enumerator is guaranteed to output canonical
  /// configurations
  bool canonical_guarantee() const;

  static const std::string enumerator_name;

 private:
  /// Implements increment
  void increment() override;

  // -- Unique -------------------

  /// Returns true if current() is valid for output
  ///
  /// When enumerating on all sites:
  /// - Require output configurations are primitive and canonical
  /// When enumerating on a subset of sites:
  /// - Require output configurations are primitive (but not canonical)
  ///
  bool _current_is_valid_for_output() const;

  /// Site index to enumerate on
  std::set<Index> m_site_index_selection;

  /// Counter over allowed occupation indices on sites in m_site_index_selection
  Counter<std::vector<int> > m_counter;

  /// The current configuration
  notstd::cloneable_ptr<Configuration> m_current;

  /// True if enumerating on a subset of supercell sites
  bool m_enumerate_on_a_subset_of_supercell_sites;
};

/** @}*/
}  // namespace CASM

#endif
