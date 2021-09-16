#ifndef CASM_config_UnitCellCoordRep_impl
#define CASM_config_UnitCellCoordRep_impl

#include "casm/configuration/UnitCellCoordRep.hh"
#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {
namespace config {

// --- Inline definitions ---

/// \brief Apply symmetry to UnitCellCoord
///
/// \param rep BasisPermuteRep representation of the symmetry operation
/// \param integral_site_coordinate Coordinate being transformed
///
inline UnitCellCoord &apply(BasisPermuteRep const &rep,
                            UnitCellCoord &integral_site_coordinate) {
  integral_site_coordinate = copy_apply(rep, integral_site_coordinate);
  return integral_site_coordinate;
}

/// \brief Apply symmetry to UnitCellCoord
///
/// \param rep BasisPermuteRep representation of the symmetry operation
/// \param integral_site_coordinate Coordinate being transformed
///
inline UnitCellCoord copy_apply(BasisPermuteRep const &rep,
                                UnitCellCoord const &integral_site_coordinate) {
  UnitCell unitcell_indices =
      rep.point_matrix * integral_site_coordinate.unitcell_indices +
      rep.unitcell_indices[integral_site_coordinate.sublattice_index];
  Index sublattice_index =
      rep.sublattice_index[integral_site_coordinate.sublattice_index];
  return UnitCellCoord(sublattice_index, unitcell_indices);
}

}  // namespace config
}  // namespace CASM

#endif
