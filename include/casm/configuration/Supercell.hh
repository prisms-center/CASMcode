#ifndef CASM_config_Supercell
#define CASM_config_Supercell

#include "casm/configuration/Prim.hh"
#include "casm/configuration/SupercellSymInfo.hh"
#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/crystallography/Superlattice.hh"

namespace CASM {
namespace config {

/// \brief Specifies all the structural and symmetry information common for all
/// configurations with the same supercell. All members are const.
struct Supercell {
  Supercell(std::shared_ptr<Prim const> const &_prim,
            Superlattice const &_superlattice);

  /// \brief Species the primitive crystal structure (lattice and basis) and
  /// allowed degrees of freedom (DoF), and also symmetry representations
  /// used for all configurations with the same prim
  std::shared_ptr<Prim const> const prim;

  /// Couples the primitive lattice to the supercell lattice, and knows the
  /// transformation matrix
  Superlattice const superlattice;

  /// \brief Converts between ijk (UnitCell) values and their corresponding
  /// index in an unrolled vector
  UnitCellIndexConverter const unitcell_index_converter;

  /// \brief Converts between bijk (UnitCellCoord) values and their
  /// corresponding linear index
  UnitCellCoordIndexConverter const unitcellcoord_index_converter;

  /// \brief Holds symmetry representations used for all configurations with
  /// the same supercell
  SupercellSymInfo const sym_info;
};

}  // namespace config
}  // namespace CASM

#endif
