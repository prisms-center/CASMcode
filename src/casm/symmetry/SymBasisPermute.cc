#include "casm/symmetry/SymBasisPermute.hh"

#include "casm/basis_set/DoF.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {

// ---- SymBasisPermute Definitions --------------------

/// Construct SymBasisPermute
SymBasisPermute::SymBasisPermute(SymOp const &_op, Lattice const &_lat,
                                 std::vector<UnitCellCoord> const &_ucc_permute)
    : m_ucc_permute(_ucc_permute),
      m_point_mat(lround(cart2frac(_op.matrix(), _lat))) {
  if (!is_integer(cart2frac(_op.matrix(), _lat), _lat.tol())) {
    throw std::runtime_error(
        std::string("Error in 'SymBasisPermute(const SymOp& op, const "
                    "StrucType& struc, double tol)'\n") +
        "  Could not get integer point transformation matrix.");
  }
}

}  // namespace CASM
