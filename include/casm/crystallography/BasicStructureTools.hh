#ifndef BASICSTRUCTURETOOLS_HH
#define BASICSTRUCTURETOOLS_HH

#include <vector>

#include "casm/crystallography/Lattice.hh"
#include "casm/external/Eigen/Core"
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {
namespace xtal {
class Site;
class Coordinate;
class BasicStructure;
struct SymOp;

/// return basis index of site that is same type and within distance 'tol' (in
/// Angstr) of test_site If no such site exists in basis, return the size of the
/// basis
Index find_index(const std::vector<Site> &basis, const Site &test_site,
                 double tol);

/// Returns true if the structure describes a crystal primitive cell
/// i.e., no translation smaller than a lattice vector can map the structure
/// onto itself
bool is_primitive(const BasicStructure &struc, double tol = TOL);

/// Returns the smallest possible tiling unit of the given structure
BasicStructure make_primitive(const BasicStructure &non_primitive_struc,
                              double tol = TOL);

/// Calculates the rotation angle and axis of a symmetry operation. This
/// function is almost exactly identical to the constructor of SymInfo::SymInfo
std::pair<double, Eigen::Vector3d> calc_rotation_angle_and_axis(
    const SymOp &op, const Lattice &lat);

/// Sort a factor group based on a lexicographical comparison of (-det, -trace,
/// angle, axis, tau)
void sort_factor_group(std::vector<SymOp> &factor_group, const Lattice &lat);

/// Create the factor group of the given structure. If the structure has no
/// degrees of freedom affected by time reversal, time reversal is ignored.
/// Otherwise symmetry operations are checked for time reversal
std::vector<SymOp> make_factor_group(const BasicStructure &struc,
                                     double tol = TOL);

/// Create the permtuation group of the given structure. The permutation group
/// can be used to identify how the basis of this structure are transformed
/// under the application of a symmetry operation.
std::vector<Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, Index>>
make_permutation_representation(const xtal::BasicStructure &struc,
                                const std::vector<SymOp> &factor_group);

// TODO: Rename to symmetrize_basis and have it take vector<Site> and Lattice?
// seems like a symmetrize routine that takes a structure should also symmetrize
// the lattice.
/// Given a symmetry group, the basis of the structure will have
/// each operation applied to it. The resulting set of basis
/// from performing these operations will be averaged out,
/// yielding a new average basis.
BasicStructure symmetrize(const BasicStructure &structure,
                          const std::vector<SymOp> &enforced_group);

/// Given an integer transformation matrix, create a superstructure whose
/// lattice is the product of the original lattice and the transformation
/// matrix.
template <typename IntegralType, int Options = 0>
BasicStructure make_superstructure(
    const BasicStructure &tiling_unit,
    const Eigen::Matrix<IntegralType, 3, 3, Options> &transformation_matrix);

}  // namespace xtal
}  // namespace CASM

#endif
