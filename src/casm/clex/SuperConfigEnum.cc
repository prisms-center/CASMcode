#include "casm/clex/SuperConfigEnum.hh"

#include "casm/app/casm_functions.hh"
#include "casm/app/enum.hh"
#include "casm/app/io/json_io_impl.hh"
#include "casm/clex/ConfigDoFTools.hh"
#include "casm/clex/ConfigEnumByPermutation.hh"
#include "casm/clex/FillSupercell.hh"
#include "casm/clex/FilteredConfigIterator.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/database/Selection.hh"

namespace CASM {

std::string SuperConfigEnum::name() const { return enumerator_name; }

const std::string SuperConfigEnum::enumerator_name = "SuperConfigEnum";

void SuperConfigEnum::_init() {
  // check that all sub-config have same supercell
  for (auto it = sub_config().begin(); it != sub_config().end(); ++it) {
    if (&it->supercell() != &(sub_config().begin()->supercell())) {
      throw std::runtime_error(
          "Error constructing SuperConfigEnum: "
          "Sub-Configurations with different Supercells");
    }
  }
  m_sub_scel = &(m_sub_config.begin()->supercell());

  // tiling of sub-config supercell lattice into super-config supercell lattice
  m_unitcell_index_converter =
      notstd::make_cloneable<xtal::UnitCellIndexConverter>(
          xtal::make_transformation_matrix_to_super(
              _sub_supercell().lattice(), _target_supercell().lattice(),
              _target_supercell().lattice().tol()));

  // initialize 'm_counter' to count over all possible sub-config on
  // each lattice point
  Index total_sites = this->_unitcell_index_converter().total_sites();
  m_counter =
      Counter<Array<int> >(Array<int>(total_sites, 0),
                           Array<int>(total_sites, sub_config().size() - 1),
                           Array<int>(total_sites, 1));

  // save indices for mapping sub-config dof values to super config dof values
  //
  // for occupation:
  //     m_current->occ(m_index_map[i][j]) = m_sub_scel[m_counter()[i]].occ(j)
  //
  // and similar for all other site DoF:
  //     m_current->local_dofs()[dof_key].values().col(m_index_map[i][j]) =
  //       m_sub_scel[m_counter()[i]].local_dofs()[dof_key].values().col(j)
  //
  m_index_map.resize(total_sites);
  for (int i = 0; i < total_sites; ++i) {
    UnitCell ref = _sub_supercell().transf_mat().cast<Index>() *
                   this->_unitcell_index_converter()(i);
    for (int j = 0; j < _sub_supercell().num_sites(); ++j) {
      UnitCellCoord uccord = _sub_supercell().uccoord(j) + ref;
      Index linear_index = _target_supercell().linear_index(uccord);
      m_index_map[i].push_back(linear_index);
    }
  }

  m_current = notstd::make_cloneable<Configuration>(_target_supercell());

  _initialize(&(*m_current));
  _fill(counter(), *m_current);

  // Make sure that current() satisfies requested conditions
  if (!_check_current()) {
    increment();
  }

  // set step to 0
  if (valid()) {
    _set_step(0);
  }
  m_current->set_source(source(step()));
}

// **** Mutators ****
// increment m_current and return a reference to it
void SuperConfigEnum::increment() {
  bool is_valid_config{false};

  while (!is_valid_config && ++m_counter) {
    _fill(counter(), *m_current);
    is_valid_config = _check_current();
  }

  if (m_counter.valid()) {
    _increment_step();
  } else {
    _invalidate();
  }
}

/// Returns true if current() satisfies requested conditions
bool SuperConfigEnum::_check_current() const { return true; }

/// Fill DoF from sub_config into a super configuration
///
/// \param counter_val The index of the sub_config on each lattice point
/// \param super_config The super configuration to set the DoF values
///
void SuperConfigEnum::_fill(Array<int> const &counter_val,
                            Configuration &super_config) {
  double xtal_tol = super_config.supercell().prim().lattice().tol();
  super_config.configdof() = make_configdof(super_config.supercell(), xtal_tol);

  for (Index i = 0; i < this->_unitcell_index_converter().total_sites(); ++i) {
    Configuration const &sub_config_i = _sub_config()[counter_val[i]];
    for (Index j = 0; j < _sub_supercell().num_sites(); ++j) {
      // copy site DoF
      super_config.set_occ(m_index_map[i][j], sub_config_i.occ(j));

      for (auto &pair : super_config.configdof().local_dofs()) {
        auto const &dof_key = pair.first;

        auto &local_dof = super_config.configdof().local_dof(dof_key);

        local_dof.site_value(m_index_map[i][j]) =
            sub_config_i.configdof().local_dof(dof_key).site_value(j);
      }
    }
  }
}

}  // namespace CASM
