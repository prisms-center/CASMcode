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

/// \brief Conditionally true for ConfigEnumAllOccupations
template <>
bool is_guaranteed_for_database_insert(
    ConfigEnumAllOccupations const &enumerator);

/// Enumerate over all possible occupations on particular sites in a
/// Configuration
///
class ConfigEnumAllOccupations : public InputEnumeratorBase<Configuration> {
  // -- Required members -------------------

 public:
  /// \brief Construct with a ConfigEnumInput, specifying which sites
  ///     to enumerate on and which to keep fixed
  ConfigEnumAllOccupations(ConfigEnumInput const &config_enum_input);

  /// \brief Constructor allowing direct control of whether
  ///     non-primitive and non-canonical Configuration are enumerated
  ConfigEnumAllOccupations(ConfigEnumInput const &config_enum_input,
                           bool primitive_only, bool canonical_only);

  std::string name() const override;

  /// \brief Returns true if enumerator is guaranteed to output
  ///     primitive & canonical configurations only
  bool primitive_canonical_guarantee() const;

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

  /// True if only enumerating primitive configurations
  bool m_primitive_only;

  /// True if only enumerating canonical configurations
  bool m_canonical_only;
};

/** @}*/
}  // namespace CASM

#endif
