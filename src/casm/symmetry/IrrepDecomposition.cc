#include "casm/symmetry/IrrepDecomposition.hh"

#include <iostream>

#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/symmetry/IrrepDecompositionImpl.hh"

namespace CASM {

namespace SymRepTools_v2 {

IrrepInfo::IrrepInfo(Eigen::MatrixXcd _trans_mat, Eigen::VectorXcd _characters)
    : trans_mat(std::move(_trans_mat)),
      characters(std::move(_characters)),
      complex(!almost_zero(trans_mat.imag())),
      pseudo_irrep(false),
      index(0) {}

/// Construct a "dummy" IrrepInfo with user specified transformtion matrix
///
/// The "dummy" IrrepInfo is constructed with specified transformtion matrix
/// and character vector of [(dim,0)] where 'dim' is the dimension of irrep
/// (number of rows of `trans_mat`)
IrrepInfo make_dummy_irrep_info(Eigen::MatrixXcd const &trans_mat) {
  Eigen::VectorXcd tchar(1);
  tchar(0) = std::complex<double>(double(trans_mat.rows()), 0.);
  return IrrepInfo(trans_mat, tchar);
}

/// Construct a "dummy" IrrepInfo with user specified transformtion matrix
///
/// The "dummy" IrrepInfo is constructed with specified transformtion matrix
/// and character vector of [(dim,0)] where 'dim' is the dimension of irrep
/// (number of rows of `trans_mat`)
IrrepInfo make_dummy_irrep_info(Eigen::MatrixXd const &trans_mat) {
  Eigen::VectorXcd tchar(1);
  tchar(0) = std::complex<double>(double(trans_mat.rows()), 0.);
  return IrrepInfo(trans_mat.template cast<std::complex<double>>(), tchar);
}

/// \brief Assumes that irreps are real, and concatenates their individual
/// trans_mats to form larger trans_mat
Eigen::MatrixXd full_trans_mat(std::vector<IrrepInfo> const &irreps) {
  Index row = 0;
  Index col = 0;
  for (auto const &irrep : irreps) {
    col = irrep.vector_dim();
    row += irrep.irrep_dim();
  }
  Eigen::MatrixXd trans_mat(row, col);
  row = 0;
  for (auto const &irrep : irreps) {
    trans_mat.block(row, 0, irrep.irrep_dim(), irrep.vector_dim()) =
        irrep.trans_mat.real();
    row += irrep.irrep_dim();
  }
  return trans_mat;
}

/// IrrepDecomposition constructor
///
/// \param rep Full space matrix representation (rep[0].rows() ==
///     init_subspace.rows())
/// \param head_group Group for which the irreps are to be found
/// \param init_subspace Input subspace in which irreps are to be found. Will be
///     expanded (column space increased) by application of rep and
///     orthogonalization to form an invariant subspace (i.e. column space
///     dimension is not increased by application of elements in head_group)
/// \param _cyclic_subgroups Cyclic subgroups of head_group. Cyclic subgropus
///     are those formed by repeated application of a single element. Used for
///     symmetrization of the irrep subspaces.
/// \param _all_subgroups All subgroups of head_group. Used for
///     symmetrization of the irrep subspaces if symmetrization using
///     _cyclic_subgroups fails.
/// \param allow_complex If true, all irreps may be complex-valued, if false,
///     complex irreps are combined to form real representations
///
IrrepDecomposition::IrrepDecomposition(
    MatrixRep const &_fullspace_rep, GroupIndices const &_head_group,
    Eigen::MatrixXd const &init_subspace,
    GroupIndicesOrbitVector const &_cyclic_subgroups,
    GroupIndicesOrbitVector const &_all_subgroups, bool allow_complex)
    : fullspace_rep(_fullspace_rep),
      head_group(_head_group),
      cyclic_subgroups(_cyclic_subgroups),
      all_subgroups(_all_subgroups) {
  using namespace IrrepDecompositionImpl;

  Index dim = fullspace_rep[0].rows();

  // 1) Expand subspace by application of group, and orthonormalization
  subspace = make_invariant_space(fullspace_rep, head_group, init_subspace);

  // 2) Perform irrep_decomposition
  // In some cases the `irrep_decomposition` method does not find all irreps.
  // As long as it finds at least one, this loop will try again in the remaining
  // subspace.
  Eigen::MatrixXd subspace_i = subspace;
  Eigen::MatrixXd finished_subspace = make_kernel(subspace);
  while (finished_subspace.cols() != dim) {
    // Irreps are found in a subspace specified via the subspace matrix rep
    MatrixRep subspace_rep_i = make_subspace_rep(fullspace_rep, subspace_i);
    std::vector<IrrepInfo> subspace_irreps_i =
        irrep_decomposition(subspace_rep_i, head_group, allow_complex);

    // If not irreps found in the subspace, this method has failed
    // If the irreps do not span the whole subspace, we'll try again
    if (subspace_irreps_i.size() == 0) {
      std::stringstream msg;
      msg << "Error in IrrepDecomposition: failed to find all irreps";
      throw std::runtime_error(msg.str());
    }

    // Symmetrize all the irreps that were found
    std::vector<IrrepInfo> symmetrized_subspace_irreps_i =
        symmetrize_irreps(subspace_rep_i, head_group, subspace_irreps_i,
                          cyclic_subgroups, all_subgroups);
    // Transform the irreps trans_mat to act on vectors in the fullspace
    std::vector<IrrepInfo> symmetrized_fullspace_irreps_i =
        make_fullspace_irreps(symmetrized_subspace_irreps_i, subspace_i);
    // Save the new fullspace irreps
    for (auto const irrep : symmetrized_fullspace_irreps_i) {
      irreps.push_back(irrep);
    }

    // Combine the irrep spaces and add to finished_subspace
    Eigen::MatrixXd finished_subspace_i =
        full_trans_mat(symmetrized_fullspace_irreps_i).adjoint();
    finished_subspace = extend(finished_subspace, finished_subspace_i);

    // If not all irreps have been found, try again in remaining space
    subspace_i = make_kernel(finished_subspace);
  }

  // 3) Combine to form symmetry adapted subspace
  symmetry_adapted_subspace = full_trans_mat(irreps).adjoint();
}

}  // namespace SymRepTools_v2

}  // namespace CASM

// for `make_irrep_decomposition` only
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymGroupRep.hh"

namespace CASM {

/// Make an IrrepDecompotion using CASM::SymGroupRep and CASM::SymGroup
SymRepTools_v2::IrrepDecomposition make_irrep_decomposition(
    SymGroupRep const &rep, SymGroup const &head_group,
    Eigen::MatrixXd const &init_subspace, bool allow_complex) {
  SymRepTools_v2::MatrixRep matrix_rep;
  for (Index i = 0; i < rep.size(); ++i) {
    matrix_rep.push_back(*rep[i]->MatrixXd());
  }
  SymRepTools_v2::GroupIndices head_group_indices;
  for (SymOp const &op : head_group) {
    head_group_indices.insert(op.index());
  }
  SymRepTools_v2::GroupIndicesOrbitVector cyclic_subgroups =
      head_group.small_subgroups();
  SymRepTools_v2::GroupIndicesOrbitVector all_subgroups =
      head_group.subgroups();
  SymRepTools_v2::IrrepDecomposition irrep_decomposition{
      matrix_rep,       head_group_indices, init_subspace,
      cyclic_subgroups, all_subgroups,      allow_complex};
  return irrep_decomposition;
}

}  // namespace CASM
