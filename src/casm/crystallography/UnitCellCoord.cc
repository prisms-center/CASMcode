#include "casm/crystallography/UnitCellCoord.hh"

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/Superlattice.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
namespace xtal {

// TODO: Make this take a tolerance
UnitCell UnitCell::from_coordinate(Coordinate const &lattice_point) {
  Eigen::Vector3l rounded_lattice_point = lround(lattice_point.const_frac());
  Eigen::Vector3d round_error =
      lattice_point.const_frac() - rounded_lattice_point.cast<double>();

  /* auto unscaled_round_error = lattice_point.const_frac() -
   * rounded_lattice_point.cast<double>(); */
  /* auto round_error =
   * lattice_point.lattice().lat_column_mat()*unscaled_round_error; //scale the
   * error proportionally to the length of the vector */

  double tol = lattice_point.lattice().tol();
  if (!almost_zero(round_error(0), tol) || !almost_zero(round_error(1), tol) ||
      !almost_zero(round_error(2), tol)) {
    std::cerr << round_error << std::endl;
    UnitCell::_throw_large_rounding_error();
  }

  return UnitCell(rounded_lattice_point);
}

UnitCell UnitCell::from_cartesian(const Eigen::Vector3d &cartesian_coord,
                                  const Lattice &tiling_unit) {
  Coordinate as_coord(cartesian_coord, tiling_unit, CART);
  return UnitCell::from_coordinate(as_coord);
}

Coordinate UnitCell::coordinate(const Lattice &tiling_unit) const {
  return Coordinate(this->cast<double>(), tiling_unit, FRAC);
}

// TODO: Make this take a tolerance
UnitCell UnitCell::reset_tiling_unit(const Lattice &current_tiling_unit,
                                     const Lattice &new_tiling_unit) const {
  auto as_coord = this->coordinate(current_tiling_unit);
  as_coord.set_lattice(new_tiling_unit, CART);
  return UnitCell::from_coordinate(as_coord);
}

//************************************************************************************************************//

UnitCellCoord UnitCellCoord::from_coordinate(const PrimType &prim,
                                             const Coordinate &coord,
                                             double tol) {
  Coordinate coord_in_prim(prim.lattice());
  coord_in_prim.cart() = coord.cart();

  for (Index b = 0; b < prim.basis().size(); ++b) {
    if (coord_in_prim.min_dist(prim.basis()[b]) < tol) {
      UnitCell coord_unitcell(
          lround(coord_in_prim.const_frac() - prim.basis()[b].const_frac()));
      return UnitCellCoord(b, coord_unitcell);
    }
  }

  throw std::runtime_error(
      "Error constructing UnitCellCoord. No basis site could be found within "
      "the given tolerance.");
}

/// \brief Get corresponding coordinate
Coordinate UnitCellCoord::coordinate(const PrimType &prim) const {
  return site(prim);
}

bool UnitCellCoord::_is_compatible_with_prim(const PrimType &prim) const {
  return this->sublattice() < prim.basis().size();
}

void UnitCellCoord::_throw_incompatible_primitive_cell() {
  throw std::runtime_error(
      "Error in UnitCellCoord. Sublattice index out of range.");
}

/// \brief Get corresponding site
Site UnitCellCoord::site(const PrimType &prim) const {
  if (!this->_is_compatible_with_prim(prim)) {
    UnitCellCoord::_throw_incompatible_primitive_cell();
  }
  return prim.basis()[this->sublattice()] +
         this->unitcell().coordinate(prim.lattice());
}

/// \brief Get reference to corresponding sublattice site in the unit structure
const Site &UnitCellCoord::sublattice_site(const PrimType &prim) const {
  if (!this->_is_compatible_with_prim(prim)) {
    UnitCellCoord::_throw_incompatible_primitive_cell();
  }
  return prim.basis()[sublattice()];
}

Coordinate make_superlattice_coordinate(const UnitCell &ijk,
                                        const Lattice &tiling_unit,
                                        const Lattice &superlattice) {
  Coordinate tcoord = ijk.coordinate(tiling_unit);
  tcoord.set_lattice(superlattice, CART);
  return tcoord;
}

Coordinate make_superlattice_coordinate(const UnitCell &ijk,
                                        const Superlattice &superlattice) {
  return make_superlattice_coordinate(ijk, superlattice.prim_lattice(),
                                      superlattice.superlattice());
}
}  // namespace xtal

}  // namespace CASM

namespace std {
std::size_t hash<CASM::xtal::UnitCell>::operator()(
    const CASM::xtal::UnitCell &value) const {
  std::size_t seed = value.size();
  for (int i = 0; i < 3; ++i) {
    seed ^= value(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  return seed;
}

std::size_t hash<CASM::xtal::UnitCellCoord>::operator()(
    const CASM::xtal::UnitCellCoord &value) const {
  const auto &uc = value.unitcell();
  auto seed = std::hash<CASM::xtal::UnitCell>()(uc);
  seed ^= value.sublattice() + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  return seed;
}
}  // namespace std
