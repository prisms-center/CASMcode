#ifndef CASM_symmetry_IrrepDecomposition
#define CASM_symmetry_IrrepDecomposition

#include <set>

#include "casm/container/multivector.hh"
#include "casm/external/Eigen/Core"
#include "casm/global/definitions.hh"

namespace CASM {

namespace SymRepTools_v2 {

typedef std::vector<Eigen::MatrixXd> MatrixRep;
typedef std::set<Index> GroupIndices;
typedef std::set<GroupIndices> GroupIndicesOrbit;
typedef std::vector<GroupIndicesOrbit> GroupIndicesOrbitVector;

struct IrrepInfo {
  /// \brief Construct an IrrepInfo with transformation matrix and vector of
  /// irreducible characters
  IrrepInfo(Eigen::MatrixXcd _trans_mat, Eigen::VectorXcd _characters);

  /// Dimension of irreducible vector space (less than or equal to vector_dim())
  Index irrep_dim() const { return trans_mat.rows(); }

  // Dimension of initial vector space (greater than or equal to irrep_dim())
  Index vector_dim() const { return trans_mat.cols(); }

  /// irrep_dim() x vector_dim() matrix that transforms a vector from the
  /// initial vector space into a vector in the irreducible vector space
  Eigen::MatrixXcd trans_mat;

  /// vector containing complex character of each group operation's action on
  /// the irreducible vector space
  Eigen::VectorXcd characters;

  /// true if any character has non-zero imaginary component, false otherwise
  bool complex;

  /// true if irrep is real but was created as direct sum of two complex irreps
  /// in this case, the 'irrep' is reducible, but this is the most-reduced
  /// representation that can still have real basis vectors
  bool pseudo_irrep;

  /// sequentially-assigned index used to distinguish between identical irreps
  /// irreps are identical if they have the same character vectors
  Index index;

  /// Vectors in the initial vector space that correspond to high-symmetry
  /// directions in the irreducible vector space. directions[i] is the i'th
  /// orbit of equivalent high-symmetry directions and directions[i].size() is
  /// the symmetric multiplicity of a direction in that orbit
  std::vector<std::vector<Eigen::VectorXd>> directions;
};

/// Construct a "dummy" IrrepInfo with user specified transformtion matrix
IrrepInfo make_dummy_irrep_info(Eigen::MatrixXcd const &trans_mat);

/// Construct a "dummy" IrrepInfo with user specified transformtion matrix
IrrepInfo make_dummy_irrep_info(Eigen::MatrixXd const &trans_mat);

/// \brief Assumes that irreps are real, and concatenates their individual
/// trans_mats to form larger trans_mat
Eigen::MatrixXd full_trans_mat(std::vector<IrrepInfo> const &irreps);

/// Performs irreducible subspace construction and symmetrization
struct IrrepDecomposition {
  /// IrrepDecomposition constructor
  IrrepDecomposition(MatrixRep const &_fullspace_rep,
                     GroupIndices const &_head_group,
                     Eigen::MatrixXd const &_init_subspace,
                     GroupIndicesOrbitVector const &_cyclic_subgroups,
                     GroupIndicesOrbitVector const &_all_subgroups,
                     bool allow_complex);

  /// Full space matrix representation of head_group
  ///
  /// fullspace_rep[i].rows() == full space dimension
  /// fullspace_rep[i].cols() == full space dimension
  MatrixRep fullspace_rep;

  /// Group (as indices into fullspace_rep) used to find irreps
  GroupIndices head_group;

  // Cyclic subgroups of head_group (only first group in each orbit is used)
  GroupIndicesOrbitVector cyclic_subgroups;

  // All subgroups of head_group (only first group in each orbit is used)
  GroupIndicesOrbitVector all_subgroups;

  /// Space in which to find irreducible subspaces. This space is formed by
  /// expanding `init_subspace`, if necessary, by application of `rep` and
  /// orthogonalization to form an invariant subspace (i.e. column space does
  /// not change upon application of elements in head_group)
  ///
  /// subspace.rows() == full space dimension
  /// subspace.cols() == dimension of invariant subspace
  Eigen::MatrixXd subspace;

  /// Irreducible spaces, symmetrized using `make_irrep_special_directions` and
  /// `make_irrep_symmetrizer_matrix` to align the irreducible space bases along
  /// high symmetry directions. Irreps are found in the `subspace` and then
  /// converted to full space dimension (meaning `irrep[i].vector_dim() == full
  /// space dimension` and `sum_i irrep[i].irrep_dim() == subspace.cols()`).
  std::vector<IrrepInfo> irreps;

  /// Symmetry adapted subspace
  ///
  /// symmetry_adapted_subspace.rows() == full space dimension
  /// symmetry_adapted_subspace.cols() == subspace.cols()
  Eigen::MatrixXd symmetry_adapted_subspace;
};

}  // namespace SymRepTools_v2

class SymGroup;
class SymGroupRep;

/// Make an IrrepDecompotion using CASM::SymGroupRep and CASM::SymGroup
SymRepTools_v2::IrrepDecomposition make_irrep_decomposition(
    SymGroupRep const &rep, SymGroup const &head_group,
    Eigen::MatrixXd const &_init_subspace, bool allow_complex);

}  // namespace CASM

#endif
