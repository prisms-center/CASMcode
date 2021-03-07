#ifndef CASM_symmetry_IrrepDecompositionImpl
#define CASM_symmetry_IrrepDecompositionImpl

#include "casm/symmetry/IrrepDecomposition.hh"

namespace CASM {

namespace SymRepTools_v2 {

namespace IrrepDecompositionImpl {

/// \brief Round entries that are within tol of being integer to that integer
/// value
Eigen::MatrixXcd prettyc(const Eigen::MatrixXcd &M);

/// \brief Round entries that are within tol of being integer to that integer
/// value
Eigen::MatrixXd pretty(const Eigen::MatrixXd &M);

// these methods are implementation details of IrrepDecomposition, they are
// included here for easier testing

Eigen::MatrixXd real_I(Index rows, Index cols);

Eigen::MatrixXd real_Zero(Index rows, Index cols);

Eigen::MatrixXcd complex_I(Index rows, Index cols);

Eigen::MatrixXcd complex_Zero(Index rows, Index cols);

/// Counts over pairs of columns (i,j), where j>=i, and phase=[1, i]
/// - skips i==j when phase==i
struct CommuterParamsCounter {
  CommuterParamsCounter();

  void reset(Eigen::MatrixXcd const &kernel);

  bool valid() const;

  bool increment();

  Index kernel_column_i;
  Index kernel_column_j;
  std::complex<double> phase;

 private:
  bool m_valid;
  Index m_max_cols;
  Index m_phase_index;  // 0: 1+0i,  1: 0+i
};

Eigen::MatrixXcd make_commuter(CommuterParamsCounter const &params,
                               MatrixRep const &rep,
                               GroupIndices const &head_group,
                               Eigen::MatrixXcd const &kernel);

Eigen::MatrixXcd make_kernel(Eigen::MatrixXcd const &subspace);
Eigen::MatrixXd make_kernel(Eigen::MatrixXd const &subspace);

Index find_end_of_equal_eigenvalues(Index begin,
                                    Eigen::VectorXd const &eigenvalues);

/// Make an irreducible subspace
Eigen::MatrixXcd make_irrep_subspace(Eigen::MatrixXcd const &KV_matrix,
                                     Index begin, Index end,
                                     bool allow_complex);

/// Calculate character for all matrices in rep
Eigen::VectorXcd make_characters(std::vector<Eigen::MatrixXcd> const &rep);

/// Calculate character for all matrices in rep
Eigen::VectorXd make_characters(std::vector<Eigen::MatrixXd> const &rep);

/// Check if approximately zero outside block along diagonal
bool make_is_block_diagonal(std::vector<Eigen::MatrixXcd> const &rep,
                            Index begin, Index end, double tol);

/// Find characters for block in range [begin, end)
Eigen::VectorXcd make_irrep_characters(std::vector<Eigen::MatrixXcd> const &rep,
                                       Index begin, Index end);

double make_squared_norm(Eigen::VectorXcd const &characters);

double make_squared_norm(Eigen::VectorXd const &characters);

std::complex<double> frobenius_product(Eigen::MatrixXcd const &matrix);

Eigen::MatrixXcd normalize_commuter(Eigen::MatrixXcd const &commuter);

/// Return true if space_A is extended by space_B
bool is_extended_by(Eigen::MatrixXcd const &space_A,
                    Eigen::MatrixXcd const &space_B);

/// Return matrix combining columns of space_A and space_B
Eigen::MatrixXcd extend(Eigen::MatrixXcd const &space_A,
                        Eigen::MatrixXcd const &space_B);

/// Return matrix combining columns of space_A and space_B
Eigen::MatrixXd extend(Eigen::MatrixXd const &space_A,
                       Eigen::MatrixXd const &space_B);

/// Data structure used for storing and checking possible irreps
struct PossibleIrrep {
  /// \param begin: col index in `eigenvalues` and `eigenvectors` corresponding
  ///     to the beginning of this possible irreducible space
  /// \param begin, end: col indices in `eigenvalues` and `eigenvectors`
  ///     corresponding to a range of equal eigenvalues
  /// \param eigenvalues Eigenvalues of (K.adjoint() * M_new * K)
  /// \param KV_matrix K * V, where V is the eigenvector matrix of
  ///     (K.adjoint() * M_new * K)
  PossibleIrrep(Eigen::VectorXd const &eigenvalues,
                Eigen::MatrixXcd const &KV_matrix,
                std::vector<Eigen::MatrixXcd> const &transformed_rep,
                Index _head_group_size, double _tol, bool allow_complex,
                Index _begin, Index _end);

  Index head_group_size;
  double tol;
  Index begin;      /// col index in eigenvalues/eigvectors for this irrep
  Index end;        /// col index in eigenvalues/eigvectors of next irrep
  Index irrep_dim;  /// irrep dimension / number of equal eigenvalues
  bool is_block_diagonal;
  Eigen::VectorXcd characters;
  double characters_squared_norm;

  /// is_block_diagonal && characters_squared_norm ~= head_group_size;
  bool is_irrep;

  /// (K * V).block(0, begin, K.rows(), irrep_dim)
  Eigen::MatrixXcd subspace;

  /// Check if Irrep is identity
  bool is_identity() const;

  /// Check if Irrep is gerade
  bool is_gerade() const;

  bool operator<(PossibleIrrep const &other) const;
};

/// Given a new commuter matrix,
/// perform eigenvalue decomposition,
/// and construct possible irreps
std::vector<PossibleIrrep> make_possible_irreps(
    Eigen::MatrixXcd const &commuter, Eigen::MatrixXcd const &kernel,
    MatrixRep const &rep, GroupIndices const &head_group, bool allow_complex);

/// Make a vector of IrrepInfo from PossibleIrreps
std::vector<IrrepInfo> make_irrep_info(std::set<PossibleIrrep> const &irreps);

/// \brief Transforms IrrepInfo constructed for a subspace to be IrrepInfo
/// appropriate for the full space (full space dimension == subspace.rows())
///
/// \param irrep IrrepInfo constructed for a subspace (irrep.trans_mat shape
///     is (subspace.cols() x subspace.rows())
/// \param subspace The subspace that subspace_irrep was constructed for
///
/// \result IrrepInfo constructed for the full space (result.trans_mat shape
///     is (subspace.cols() x subspace.cols())
///
IrrepInfo subspace_to_full_space(IrrepInfo const &subspace_irrep,
                                 Eigen::MatrixXd const &subspace);

/// Construct a "dummy" IrrepInfo with user specified transformtion matrix
///
/// The "dummy" IrrepInfo is constructed with specified transformtion matrix
/// and character vector of [(dim,0)] where 'dim' is the dimension of irrep
/// (number of rows of `trans_mat`)
IrrepInfo make_dummy_irrep_info(Eigen::MatrixXcd const &trans_mat);

/// Check if a representation is irreducible
///
/// A representation is irreducible if the squared norm of the characters
/// equals the group size
bool is_irrep(MatrixRep const &rep, GroupIndices const &head_group);

/// Finds irreducible subspaces that comprise an underlying subspace
std::vector<IrrepInfo> irrep_decomposition(MatrixRep const &rep,
                                           GroupIndices const &head_group,
                                           bool allow_complex);

/// Convert irreps generated for a subspace to full space dimension
std::vector<IrrepInfo> make_fullspace_irreps(
    std::vector<IrrepInfo> const &subspace_irreps,
    Eigen::MatrixXd const &subspace);

/// Expand subspace by application of group, and orthogonalize
Eigen::MatrixXd make_invariant_space(MatrixRep const &rep,
                                     GroupIndices const &head_group,
                                     Eigen::MatrixXd const &subspace);

// Create `subspace_rep`, a transformed copy of `fullspace_rep` that acts
// on coordinates with `subspace` columns as a basis. Matrices in
// `subspace_rep` are shape (subspace.cols() x subspace.cols())
MatrixRep make_subspace_rep(MatrixRep const &fullspace_rep,
                            Eigen::MatrixXd const &subspace);

/// \brief Symmetrize IrrepInfo, by finding high symmetry directions and
/// aligning the irrep subspace basis with those directions
std::vector<IrrepInfo> symmetrize_irreps(
    MatrixRep const &subspace_rep, GroupIndices const &head_group,
    std::vector<IrrepInfo> const &irreps,
    GroupIndicesOrbitVector const &cyclic_subgroups,
    GroupIndicesOrbitVector const &all_subgroups);

}  // namespace IrrepDecompositionImpl

}  // namespace SymRepTools_v2

}  // namespace CASM

#endif
