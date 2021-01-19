#ifndef CASM_SuperConfigEnum
#define CASM_SuperConfigEnum

#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/container/Counter.hh"
#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/enumerator/InputEnumerator.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

/// \brief Given a selection of Configurations, enumerate all tilings into
/// some supercell
///
/// - 'sub_config' is used to refer to the input 'building block' Configuration
///   - sub_config must all have the same Supercell
/// - 'target_scel' is used to refer to the Supercell that is being tiled by
///   sub_config
/// - enumerates all Configurations, including non-primitive and non-canonical
///
/// Notes:
/// - Make use of `is_valid_sub_configuration` and
/// `make_all_super_configurations` when preparing
///   valid constructor inputs.
///
/// \ingroup ConfigEnumGroup
///
class SuperConfigEnum : public InputEnumeratorBase<Configuration> {
  // -- Required members -------------------

 public:
  /// \brief Constructor, using all Supercell permutations
  template <typename ConfigIterator>
  SuperConfigEnum(Supercell const &_target_scel,
                  ConfigIterator sub_config_begin,
                  ConfigIterator sub_config_end);

  std::string name() const override;

  static const std::string enumerator_name;

  // -- Unique -------------------

  /// Access the sub-configurations
  std::vector<Configuration> const &sub_config() const { return m_sub_config; }

  /// \brief Access the current value of the counter
  ///
  /// - The counter indicates how the sub-configurations tile into the
  ///   super-configuration
  /// - sub_config()[counter()[i]] is tiled into the i-th lattice point location
  Array<int> const &counter() const { return m_counter(); }

 private:
  // -- Required members -------------------

  /// Implement increment
  void increment() override;

  // -- Unique -------------------

  /// Used by constructor
  void _init();

  /// Returns true if current() satisifies requested conditions
  bool _check_current() const;

  /// Access the super-config supercell
  Supercell const &_target_supercell() { return m_target_scel; }

  /// Access the sub-config supercell
  Supercell const &_sub_supercell() const { return *m_sub_scel; }

  /// Access the sub-configurations
  std::vector<Configuration> &_sub_config() { return m_sub_config; }

  /// Access the UnitCellIndexConverter
  xtal::UnitCellIndexConverter &_unitcell_index_converter() {
    return *m_unitcell_index_converter;
  }

  /// Fill DoF from sub_config into a Configuration
  ///
  /// \param summary The index of the sub_config on each lattice site
  void _fill(Array<int> const &counter_val, Configuration &config);

  /// The supercell being filled
  Supercell const &m_target_scel;

  /// A vector containing each possible sub_config
  std::vector<Configuration> m_sub_config;

  // All sub_config use the same supercell
  Supercell const *m_sub_scel;

  /// The 'current' Configuration
  notstd::cloneable_ptr<Configuration> m_current;

  Counter<Array<int> > m_counter;

  /// An ordered set of lattice points (UnitCell) that tile the sub_scel into
  /// the ref_scel
  notstd::cloneable_ptr<xtal::UnitCellIndexConverter>
      m_unitcell_index_converter;

  /// m_current->occ(m_index_map[i][j]) = m_sub_scel[i].occ(j)
  std::vector<std::vector<Index> > m_index_map;
};

/// Constructor
///
template <typename ConfigIterator>
SuperConfigEnum::SuperConfigEnum(Supercell const &_target_scel,
                                 ConfigIterator sub_config_begin,
                                 ConfigIterator sub_config_end)
    : m_target_scel(_target_scel),
      m_sub_config(sub_config_begin, sub_config_end) {
  _init();
}

}  // namespace CASM

#endif
