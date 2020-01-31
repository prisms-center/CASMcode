#ifndef CASM_SuperConfigEnum
#define CASM_SuperConfigEnum

#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/container/Counter.hh"
#include "casm/enumerator/InputEnumerator.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_SuperConfigEnum_interface();
}

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
  /// \ingroup ConfigEnumGroup
  ///
  class SuperConfigEnum : public InputEnumeratorBase<Configuration> {

    // -- Required members -------------------

  public:

    /// \brief Constructor, using all Supercell permutations
    template<typename ConfigIterator>
    SuperConfigEnum(const Supercell &_target_scel,
                    ConfigIterator _begin,
                    ConfigIterator _end);

    /// \brief Constructor
    template<typename ConfigIterator>
    SuperConfigEnum(const Supercell &_target_scel,
                    ConfigIterator _begin,
                    ConfigIterator _end,
                    PermuteIterator _perm_begin,
                    PermuteIterator _perm_end);

    std::string name() const override {
      return enumerator_name;
    }

    static const std::string enumerator_name;
    static std::string interface_help();
    static int run(const PrimClex &primclex, const jsonParser &kwargs, const Completer::EnumOption &enum_opt, EnumeratorMap const *interface_map);


    // -- Unique -------------------

    /// Access the sub-configurations
    const std::vector<Configuration> &sub_config() const {
      return m_sub_config;
    }

    /// \brief Access the current value of the counter
    ///
    /// - The counter indicates how the sub-configurations tile into the
    ///   super-configuration
    /// - sub_config()[counter()[i]] is tiled into the i-th lattice point location
    const Array<int> &counter() const {
      return m_counter();
    }

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
    const Supercell &_target_supercell() {
      return m_target_scel;
    }

    /// Access the sub-config supercell
    const Supercell &_sub_supercell() const {
      return *m_sub_scel;
    }

    /// Access the sub-configurations
    std::vector<Configuration> &_sub_config() {
      return m_sub_config;
    }

    /// Access the UnitCellIndexConverter
    xtal::UnitCellIndexConverter &_unitcell_index_converter() {
      return *m_unitcell_index_converter;
    }

    /// Fill DoF from sub_config into a Configuration
    ///
    /// \param summary The index of the sub_config on each lattice site
    void _fill(const Array<int> &counter_val, Configuration &config);

    const PermuteIterator &_perm_begin() const {
      return m_perm_begin;
    }
    const PermuteIterator &_perm_end() const {
      return m_perm_end;
    }


    /// The supercell being filled
    const Supercell &m_target_scel;

    /// A vector containing each possible sub_config
    std::vector<Configuration> m_sub_config;

    PermuteIterator m_perm_begin, m_perm_end;

    // All sub_config use the same supercell
    const Supercell *m_sub_scel;

    /// The 'current' Configuration
    notstd::cloneable_ptr<Configuration> m_current;

    Counter<Array<int> > m_counter;

    /// An ordered set of lattice points (UnitCell) that tile the sub_scel into the ref_scel
    notstd::cloneable_ptr<xtal::UnitCellIndexConverter> m_unitcell_index_converter;

    /// m_current->occ(m_index_map[i][j]) = m_sub_scel[i].occ(j)
    std::vector<std::vector<Index> > m_index_map;

    bool m_has_occ;

  };

  /// \brief Constructor, using all Supercell permutations
  ///
  template<typename ConfigIterator>
  SuperConfigEnum::SuperConfigEnum(const Supercell &_target_scel,
                                   ConfigIterator _begin,
                                   ConfigIterator _end) :
    SuperConfigEnum(_target_scel,
                    _begin,
                    _end,
                    _target_scel.sym_info().permute_begin(),
                    _target_scel.sym_info().permute_end()) {}

  /// \brief Constructor
  ///
  template<typename ConfigIterator>
  SuperConfigEnum::SuperConfigEnum(const Supercell &_target_scel,
                                   ConfigIterator _begin,
                                   ConfigIterator _end,
                                   PermuteIterator _perm_begin,
                                   PermuteIterator _perm_end) :
    m_target_scel(_target_scel),
    m_sub_config(_begin, _end),
    m_perm_begin(_perm_begin),
    m_perm_end(_perm_end) {

    _init();

  }

}

#endif
