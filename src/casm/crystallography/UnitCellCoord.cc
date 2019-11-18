#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/Site.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/symmetry/SymBasisPermute.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {
  namespace xtal {

    UnitCell UnitCell::from_coordinate(Coordinate const &lattice_point) {
      return UnitCell(lround(lattice_point.const_frac()));
    }


    UnitCellCoord UnitCellCoord::from_coordinate(const PrimType &prim, const Coordinate &coord, double tol) {
      Coordinate coord_in_prim(prim.lattice());
      coord_in_prim.cart() = coord.cart();

      for(Index b = 0; b < prim.basis().size(); ++b) {
        auto coord_distance_to_basis_site = coord_in_prim - prim.basis()[b];
        auto rounded_distance = coord_distance_to_basis_site;
        rounded_distance.frac() = round(coord_distance_to_basis_site.const_frac());

        if((coord_distance_to_basis_site - rounded_distance).const_cart().norm() < tol) {
          return UnitCellCoord(b, lround(coord_distance_to_basis_site.const_frac()));
        }
      }

      throw std::runtime_error("Error constructing UnitCellCoord. No basis site could be found within the given tolerance.");
    }

    /// \brief Get corresponding coordinate
    Coordinate UnitCellCoord::coordinate(const PrimType &prim) const {
      return site(prim);
    }

    bool UnitCellCoord::_is_compatible_with_prim(const PrimType &prim) const {
      return this->sublattice() < prim.basis().size();
    }

    void UnitCellCoord::_throw_incompatible_primitive_cell() {
      throw std::runtime_error("Error in UnitCellCoord. Sublattice index out of range.");
    }

    /// \brief Get corresponding site
    Site UnitCellCoord::site(const PrimType &prim) const {
      if(!this->_is_compatible_with_prim(prim)) {
        UnitCellCoord::_throw_incompatible_primitive_cell();
      }
      return prim.basis()[sublattice()] + Coordinate(unitcell().cast<double>(), prim.lattice(), FRAC);
    }

    /// \brief Get reference to corresponding sublattice site in the unit structure
    const Site &UnitCellCoord::sublattice_site(const PrimType &prim) const {
      if(!this->_is_compatible_with_prim(prim)) {
        UnitCellCoord::_throw_incompatible_primitive_cell();
      }
      return prim.basis()[sublattice()];
    }

  } // namespace xtal

  /// \brief Print to json as [b, i, j, k]
  jsonParser &to_json(const xtal::UnitCellCoord &ucc_val, jsonParser &fill_json) {
    fill_json.put_array();
    fill_json.push_back(ucc_val.sublattice());
    fill_json.push_back(ucc_val.unitcell()(0));
    fill_json.push_back(ucc_val.unitcell()(1));
    fill_json.push_back(ucc_val.unitcell()(2));

    return fill_json;
  }

  /// \brief Read from json [b, i, j, k], assuming fill_value.unit() is already set
  void from_json(UnitCellCoord &fill_value, const jsonParser &read_json) {

    auto b = read_json[0].get<Index>();
    auto i = read_json[1].get<Index>();
    auto j = read_json[2].get<Index>();
    auto k = read_json[3].get<Index>();

    fill_value = UnitCellCoord(b, i, j, k);

    return;
  }
} // namespace CASM
