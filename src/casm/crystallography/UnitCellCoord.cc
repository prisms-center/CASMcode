#include "casm/casm_io/json/jsonParser.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
  namespace xtal {

    UnitCell UnitCell::from_coordinate(Coordinate const &lattice_point) {
      return UnitCell(lround(lattice_point.const_frac()));
    }


    UnitCellCoord UnitCellCoord::from_coordinate(const PrimType &prim, const Coordinate &coord, double tol) {
      Coordinate coord_in_prim(prim.lattice());
      coord_in_prim.cart() = coord.cart();

      for(Index b = 0; b < prim.basis().size(); ++b) {
        if(coord_in_prim.min_dist(prim.basis(b)) < tol) {
          UnitCell coord_unitcell(lround(coord_in_prim.const_frac() - prim.basis(b).const_frac()));
          return UnitCellCoord(b, coord_unitcell);
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
  void from_json(xtal::UnitCellCoord &fill_value, const jsonParser &read_json) {

    auto b = read_json[0].get<Index>();
    auto i = read_json[1].get<Index>();
    auto j = read_json[2].get<Index>();
    auto k = read_json[3].get<Index>();

    fill_value = xtal::UnitCellCoord(b, i, j, k);

    return;
  }
} // namespace CASM
