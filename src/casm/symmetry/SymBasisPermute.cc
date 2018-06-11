#include "casm/symmetry/SymBasisPermute.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/container/LinearAlgebra.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/basis_set/DoF.hh"

namespace CASM {

  // ---- SymBasisPermute Definitions --------------------

  /// Construct SymBasisPermute
  template<typename StrucType>
  SymBasisPermute::SymBasisPermute(const SymOp &op, const StrucType &struc, double tol) {
    SymOp::matrix_type frac_op(cart2frac(op.matrix(), struc.lattice()));
    if(!is_integer(frac_op, tol)) {
      throw std::runtime_error(
        std::string("Error in 'SymBasisPermute(const SymOp& op, const StrucType& struc, double tol)'\n") +
        "  Could not get integer point transformation matrix.");
    }

    m_point_mat = lround(frac_op);

    // Determine how basis sites transform from the origin unit cell
    for(int b = 0; b < struc.basis().size(); b++) {
      m_ucc_permute.push_back(UnitCellCoord(struc, CASM::copy_apply(op, struc.basis()[b]), tol));
    }
  }

  template SymBasisPermute::SymBasisPermute(const SymOp &op, const Structure &struc, double tol);
}

