#ifndef CASM_config_UnitCellCoordRep
#define CASM_config_UnitCellCoordRep

#include "casm/configuration/definitions.hh"

namespace CASM {
namespace config {

/// \brief Data structure used to transform integral site coordinates
///
/// - A symmetry operation {R,T} transforms Cartesian coordinate (x) like:
///     x' = R*x + T, where R is a point transformation matrix, and T a
///     translation vector
/// - BasisPermuteRep transforms integer fractional coordinate (u) without
///   basis like:
///    Lu' = R*L*u + T
///     u' = L.inv*R*L*u + L.inv*T, where L is the lattice, as a column vector
///     matrix
/// - For transforming basis sites, a lookup table is stored that maps
///     UnitCellCoord(struc, b, UnitCell(0,0,0)) -> UnitCellCoord'
///   which is used to set the sublat and is added to u' along with L.inv*T
/// - L.inv*R*L is stored as BasisPermuteRep::point_matrix
/// - The sublattice index after transformation is stored in
///   SymBasisPermute::sublattice_index.
/// - So an UnitCellCoord, c, is transformed by an
///   UnitCellCoordRep, rep, like:
///   \code
///   c_after.unitcell_indices =
///       rep.point_matrix * c_before.unitcell_indices +
///       rep.unitcell_indices[c_before.sublattice_index];
///   c_after.sublattice_index =
///       rep.sublattice_index[c_before.sublattice_index];
/// \endcode
/// - Generally UnitCellCoordRep are constructed for prim factor group
///   operations and additional translations, such as those in supercell factor
///   group operations, may also be necessary.
///
struct UnitCellCoordRep {
  /// \brief Sublattice index after transformation (for each sublattice)
  ///
  /// final_sublattice_index = sublattice_index[initial_sublattice_index]
  std::vector<Index> sublattice_index;

  /// \brief Unit cell indices of sites in the origin unit cell after
  /// transformation (for each sublattice)
  ///
  /// u' = point_matrix * u + unitcell_indices[initial_sublattice_index]
  std::vector<Eigen::Vector3l> unitcell_indices;

  /// \brief Integer transformation matrix for integer fractional coordinates
  ///
  /// u' = point_matrix * u + unitcell_indices[initial_sublattice_index]
  Eigen::Matrix3l point_matrix;
};

/// \brief Apply symmetry to UnitCellCoord
UnitCellCoord &apply(UnitCellCoordRep const &rep,
                     UnitCellCoord &integral_site_coordinate);

/// \brief Apply symmetry to UnitCellCoord
UnitCellCoord copy_apply(UnitCellCoordRep const &rep,
                         UnitCellCoord integral_site_coordinate);

}  // namespace config
}  // namespace CASM

#endif
