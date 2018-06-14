#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/symmetry/SymBasisPermute.hh"
#include "casm/container/LinearAlgebra.hh"
#include "casm/casm_io/jsonParser.hh"

namespace CASM {

  /// Convert lattice point a unitcell
  UnitCell make_unitcell(Coordinate const &lattice_point) {
    return UnitCell(lround(lattice_point.const_frac()));
  }

  UnitCellCoord::UnitCellCoord(const UnitType &unit, const Coordinate &coord, double tol) :
    m_unit(&unit) {
    Coordinate coord_in_unit(unit.lattice());
    coord_in_unit.cart() = coord.cart();
    for(Index b = 0; b < unit.basis().size(); ++b) {
      //Standard debugging statements when things go wrong - Please leave in
      //std::cout << "Coord" << coord_in_unit.const_frac() <<std::endl;
      //std::cout << "b" << unit.basis[b].const_frac() <<std::endl;
      auto diff = coord_in_unit - unit.basis()[b];
      //std::cout << "diff" << diff.const_frac() <<std::endl;
      Coordinate tmp = diff;
      tmp.frac() = round(diff.const_frac());

      //std::cout << "tmp" << tmp.const_frac() <<std::endl;
      //std::cout << "error is " << (diff - tmp).const_cart().norm() << "tol is " << tol<< std::endl;
      if((diff - tmp).const_cart().norm() < tol) {
        *this = UnitCellCoord(unit, b, lround(diff.const_frac()));
        return;
      }
    }

    throw std::runtime_error(
      "Error in 'UnitCellCoord(CoordType coord, const StrucType& struc, double tol)'\n"
      "  No matching basis site found.");
  }

  /// \brief Access the Lattice
  const Lattice &UnitCellCoord::lattice() const {
    return unit().lattice();
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
    if(sublat() < 0 || sublat() >= unit().basis().size()) {
      unit().print_xyz(std::cout);
      std::cerr << "CRITICAL ERROR: In BasicStructure<CoordType>::get_site(), UnitCellCoord " << *this << " is out of bounds!\n"
                << "                Cannot index basis, which contains " << unit().basis().size() << " objects.\n";
      throw std::runtime_error("Error: in 'UnitCellCoord::site()': Cannot convert UnitCellCoord to Site");
    }
    return unit().basis()[sublat()] + Coordinate(unitcell().cast<double>(), unit().lattice(), FRAC);
  }

  /// \brief Get reference to corresponding sublattice site in the unit structure
  const Site &UnitCellCoord::sublat_site() const {
    if(sublat() < 0 || sublat() >= unit().basis().size()) {
      unit().print_xyz(std::cout);
      std::cerr << "CRITICAL ERROR: In BasicStructure<CoordType>::get_site(), UnitCellCoord " << *this << " is out of bounds!\n"
                << "                Cannot index basis, which contains " << unit().basis().size() << " objects.\n";
      throw std::runtime_error("Error: in 'UnitCellCoord::site()': Cannot convert UnitCellCoord to Site");
    }
    return unit().basis()[sublat()];
  }

  UnitCellCoord &UnitCellCoord::apply_sym(const SymOp &op) {

    // transform using stored SymBasisPermute representation
    const SymBasisPermute &rep = *op.get_basis_permute_rep(unit().basis_permutation_symrep_ID());
    unitcell() = rep.matrix() * unitcell() + rep[sublat()].unitcell();
    sublat() = rep[sublat()].sublat();

    // additional translations (such as needed for supercell factor groups),
    // are stored in SymOp::integral_tau() (in cartesian coordinates)
    // this converts that to fractional coordinates and adds it to this->unitcell()
    unitcell() += lround(unit().lattice().inv_lat_column_mat() * op.integral_tau());

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
