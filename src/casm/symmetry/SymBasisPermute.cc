#include "casm/symmetry/SymBasisPermute.hh"

#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {

// ---- SymBasisPermute Definitions --------------------

/// Construct SymBasisPermute (deprecated)
///
/// Note: The method used to calculate the point operation matrix here, taking
/// into account floating point errors, is not necessary the same as used else
/// where to determine the valid symmetry operations, so it is preferrable to
/// provide `point_mat` directly.
SymBasisPermute::SymBasisPermute(SymOp const &_op, Lattice const &_lat,
                                 std::vector<UnitCellCoord> const &_ucc_permute)
    : m_ucc_permute(_ucc_permute),
      m_point_mat(lround(cart2frac(_op.matrix(), _lat))) {
  if (!is_integer(cart2frac(_op.matrix(), _lat), _lat.tol())) {
    std::cout << "fractional transformation matrix: \n"
              << cart2frac(_op.matrix(), _lat) << std::endl;
    std::cout << "m_point_mat: \n" << m_point_mat << std::endl;
    throw std::runtime_error(
        "Error constructing SymBasisPermute: Could not get integer point "
        "transformation matrix.");
  }
}

/// Construct SymBasisPermute
SymBasisPermute::SymBasisPermute(Eigen::Matrix3l const &_point_mat,
                                 std::vector<UnitCellCoord> const &_ucc_permute)
    : m_ucc_permute(_ucc_permute), m_point_mat(_point_mat) {}

}  // namespace CASM
