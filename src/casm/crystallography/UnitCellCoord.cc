#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/symmetry/SymBasisPermute.hh"


namespace CASM {

  UnitCellCoord::UnitCellCoord(const UnitType &unit, const Coordinate &coord, double tol) :
    m_unit(&unit) {
    for(Index b = 0; b < unit.basis.size(); ++b) {
      auto diff = coord - unit.basis[b];
      if(is_integer(diff.const_frac(), tol)) {
        *this = UnitCellCoord(unit, b, lround(diff.const_frac()));
        return;
      }
    }

    throw std::runtime_error(
      "Error in 'UnitCellCoord(CoordType coord, const StrucType& struc, double tol)'\n"
      "  No matching basis site found.");
  }

  /// \brief Get corresponding coordinate
  Coordinate UnitCellCoord::coordinate() const {
    return site();
  }

  UnitCellCoord::operator Coordinate() const {
    return coordinate();
  }

  /// \brief Get corresponding site
  Site UnitCellCoord::site() const {
    if(sublat() < 0 || sublat() >= unit().basis.size()) {
      unit().print_xyz(std::cout);
      std::cerr << "CRITICAL ERROR: In BasicStructure<CoordType>::get_site(), UnitCellCoord " << *this << " is out of bounds!\n"
                << "                Cannot index basis, which contains " << unit().basis.size() << " objects.\n";
      throw std::runtime_error("Error: in 'UnitCellCoord::site()': Cannot convert UnitCellCoord to Site");
    }
    return unit().basis[sublat()] + Coordinate(unitcell().cast<double>(), unit().lattice(), FRAC);
  }

  /// \brief Get reference to corresponding sublattice site in the unit structure
  const Site &UnitCellCoord::sublat_site() const {
    if(sublat() < 0 || sublat() >= unit().basis.size()) {
      unit().print_xyz(std::cout);
      std::cerr << "CRITICAL ERROR: In BasicStructure<CoordType>::get_site(), UnitCellCoord " << *this << " is out of bounds!\n"
                << "                Cannot index basis, which contains " << unit().basis.size() << " objects.\n";
      throw std::runtime_error("Error: in 'UnitCellCoord::site()': Cannot convert UnitCellCoord to Site");
    }
    return unit().basis[sublat()];
  }

  UnitCellCoord &UnitCellCoord::apply_sym(const SymOp &op) {
    const SymBasisPermute &rep = *op.get_basis_permute_rep(unit().basis_permutation_symrep_ID());
    unitcell() = rep.matrix() * unitcell() + rep[sublat()].unitcell() + (unit().lattice().inv_lat_column_mat() * op.integral_tau()).cast<long>();
    sublat() = rep[sublat()].sublat();

    return *this;
  }

  UnitCellCoord UnitCellCoord::copy_apply(const SymOp &op) const {
    UnitCellCoord result(*this);
    result.apply_sym(op);
    return result;
  }

  /// \brief Print to json as [b, i, j, k]
  jsonParser &to_json(const UnitCellCoord &ucc_val, jsonParser &fill_json) {
    fill_json.put_array();
    fill_json.push_back(ucc_val.sublat());
    fill_json.push_back(ucc_val.unitcell()(0));
    fill_json.push_back(ucc_val.unitcell()(1));
    fill_json.push_back(ucc_val.unitcell()(2));

    return fill_json;
  }

  /// \brief Read from json [b, i, j, k], using 'unit' for UnitCellCoord::unit()
  UnitCellCoord jsonConstructor<UnitCellCoord>::from_json(const jsonParser &json, const Structure &unit) {
    UnitCellCoord coord(unit);
    CASM::from_json(coord, json);
    return coord;
  }

  /// \brief Read from json [b, i, j, k], assuming fill_value.unit() is already set
  void from_json(UnitCellCoord &fill_value, const jsonParser &read_json) {

    fill_value.sublat() = read_json[0].get<Index>();
    fill_value.unitcell()(0) = read_json[1].get<Index>();
    fill_value.unitcell()(1) = read_json[2].get<Index>();
    fill_value.unitcell()(2) = read_json[3].get<Index>();

    return;
  }
}

