#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/symmetry/SymBasisPermute.hh"


namespace CASM {

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

    unitcell() = rep.matrix() * unitcell() + rep[sublat()].unitcell() + (unit().lattice().inv_lat_column_mat() * op.tau()).cast<long>();
    sublat() = rep[sublat()].sublat();

    return *this;
  }

  UnitCellCoord UnitCellCoord::copy_apply(const SymOp &op) const {
    UnitCellCoord result(*this);
    result.apply_sym(op);
    return result;
  }
}

