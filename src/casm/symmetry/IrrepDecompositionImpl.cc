#include "casm/symmetry/IrrepDecompositionImpl.hh"

#include <iostream>

#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/symmetry/Symmetrizer.hh"
#include "casm/symmetry/to_real.hh"

namespace CASM {

namespace SymRepTools_v2 {

namespace IrrepDecompositionImpl {

/// \brief Round entries that are within tol of being integer to that integer
/// value
Eigen::MatrixXcd prettyc(const Eigen::MatrixXcd &M) {
  double tol = 1e-10;
  Eigen::MatrixXcd Mp(M);
  for (int i = 0; i < M.rows(); i++) {
    for (int j = 0; j < M.cols(); j++) {
      if (std::abs(std::round(M(i, j).real()) - M(i, j).real()) < tol) {
        Mp(i, j).real(std::round(M(i, j).real()));
      }
      if (std::abs(std::round(M(i, j).imag()) - M(i, j).imag()) < tol) {
        Mp(i, j).imag(std::round(M(i, j).imag()));
      }
    }
  }
  return Mp;
}

/// \brief Round entries that are within tol of being integer to that integer
/// value
Eigen::MatrixXd pretty(const Eigen::MatrixXd &M) {
  double tol = 1e-10;
  Eigen::MatrixXd Mp(M);
  for (int i = 0; i < M.rows(); i++) {
    for (int j = 0; j < M.cols(); j++) {
      if (std::abs(std::round(M(i, j)) - M(i, j)) < tol) {
        Mp(i, j) = std::round(M(i, j));
      }
    }
  }
  return Mp;
}

// note: there are a number of floating point comparisons, currently all set
// to use CASM::TOL as the tolerance. If necessary, they could be tuned here or
// via function parameters. To find them search for TOL or "almost".

Eigen::MatrixXd real_I(Index rows, Index cols) {
  return Eigen::MatrixXd::Identity(rows, cols);
}

Eigen::MatrixXd real_Zero(Index rows, Index cols) {
  return Eigen::MatrixXd::Zero(rows, cols);
}

Eigen::MatrixXcd complex_I(Index rows, Index cols) {
  return Eigen::MatrixXcd::Identity(rows, cols);
}

Eigen::MatrixXcd complex_Zero(Index rows, Index cols) {
  return Eigen::MatrixXcd::Zero(rows, cols);
}

CommuterParamsCounter::CommuterParamsCounter() : m_valid(false) {}

void CommuterParamsCounter::reset(Eigen::MatrixXcd const &kernel) {
  m_valid = true;
  m_max_cols = kernel.cols();
  m_phase_index = 0;

  kernel_column_i = 0;
  kernel_column_j = m_max_cols - 1;
  phase = std::complex<double>(1., 0.);
}

bool CommuterParamsCounter::valid() const { return m_valid; }

bool CommuterParamsCounter::increment() {
  if (!m_valid) {
    return false;
  }

  // inner loop is column j
  --kernel_column_j;

  // next inner loop is column i
  if (kernel_column_j < kernel_column_i) {
    ++kernel_column_i;
    kernel_column_j = m_max_cols - 1;
  }

  // outer loop is phase: 1, i
  if (kernel_column_i == m_max_cols) {
    ++m_phase_index;
    kernel_column_i = 0;
    kernel_column_j = m_max_cols - 1;

    // once m_phase_index gets to 2, we've finished all possibilities
    if (m_phase_index == 2) {
      m_valid = false;
      return false;
    } else if (m_phase_index == 1) {
      phase = std::complex<double>(0., 1.);
    }
  }

  // skip i==j when phase==i
  if (kernel_column_i == kernel_column_j && m_phase_index == 1) {
    return increment();
  }

  return true;
}

Eigen::MatrixXcd make_commuter(CommuterParamsCounter const &params,
                               MatrixRep const &rep,
                               GroupIndices const &head_group,
                               Eigen::MatrixXcd const &kernel) {
  Index dim = rep[0].rows();

  // construct outer product space of kernel columns i and j
  // commuters are constructed to be self-adjoint,
  // which assures eigenvalues are real
  auto const &col_i = kernel.col(params.kernel_column_i);
  auto const &col_j = kernel.col(params.kernel_column_j);
  auto const &phase = params.phase;
  Eigen::MatrixXcd M_init = phase * col_i * col_j.adjoint() +
                            std::conj(phase) * col_j * col_i.adjoint();
  Eigen::MatrixXcd M = complex_Zero(dim, dim);

  // Reynolds operation to symmetrize:
  for (Index element_index : head_group) {
    M += rep[element_index] * M_init * rep[element_index].transpose();
  }
  return M;
}

Eigen::MatrixXcd make_kernel(Eigen::MatrixXcd const &subspace) {
  Eigen::HouseholderQR<Eigen::MatrixXcd> qr;
  qr.compute(subspace);
  return Eigen::MatrixXcd(qr.householderQ())
      .rightCols(subspace.rows() - subspace.cols());
}
Eigen::MatrixXd make_kernel(Eigen::MatrixXd const &subspace) {
  Eigen::HouseholderQR<Eigen::MatrixXd> qr;
  qr.compute(subspace);
  return Eigen::MatrixXd(qr.householderQ())
      .rightCols(subspace.rows() - subspace.cols());
}

Index find_end_of_equal_eigenvalues(Index begin,
                                    Eigen::VectorXd const &eigenvalues) {
  Index end = begin + 1;
  while (end < eigenvalues.size() &&
         almost_equal(eigenvalues(begin), eigenvalues(end), TOL)) {
    end++;
  }
  return end;
}

/// Make an irreducible subspace
///
/// Makes an orthogonalized irreducible subspace, from (K * V), where
/// - K is the kernel matrix,
/// - V is the eigenvector matrix of ( K.adjoint() * M_new * K )
/// - M_new is the new non-zero commuter matrix
///
/// \param KV_matrix Matrix (K * V)
/// \param begin, end Range of equal eigenvalues
/// \param allow_complex If true, allow subspace with complex basis vectors. If
///     false, will make a pseudo irrep subspace that combines two complex
///     irreps. In this case the irrep is reducible, but this is the most-
///     reduced representation that has real basis vectors.
///
Eigen::MatrixXcd make_irrep_subspace(Eigen::MatrixXcd const &KV_matrix,
                                     Index begin, Index end,
                                     bool allow_complex) {
  Index dim = KV_matrix.rows();

  // get columns of (K * V) corresponding to equal eigenvalues
  Eigen::MatrixXcd X = (KV_matrix).block(0, begin, dim, end - begin);
  Eigen::MatrixXcd subspace_init;
  if (allow_complex) {
    subspace_init = X;
  } else {
    subspace_init = Eigen::MatrixXcd::Zero(dim, 2 * X.cols());
    subspace_init.leftCols(X.cols()) =
        sqrt(2.0) * X.real().cast<std::complex<double>>();
    subspace_init.rightCols(X.cols()) =
        sqrt(2.0) * X.imag().cast<std::complex<double>>();
  }

  // QR decomposition

  // "it seems stupid", but Eigen::HouseholderQR is not rank revealing, and
  // Eigen::ColPivHouseholderQR permutes columns of Q
  Eigen::HouseholderQR<Eigen::MatrixXcd> qr;
  Eigen::ColPivHouseholderQR<Eigen::MatrixXcd> colqr;
  colqr.setThreshold(TOL);

  qr.compute(subspace_init);
  colqr.compute(subspace_init);
  Eigen::MatrixXcd Q = qr.householderQ();
  Eigen::MatrixXcd irrep_subspace = Q.leftCols(colqr.rank());
  return irrep_subspace;
}

/// Calculate character for all matrices in rep
Eigen::VectorXcd make_characters(std::vector<Eigen::MatrixXcd> const &rep) {
  Eigen::VectorXcd characters(rep.size());

  Index element_index = 0;
  for (Eigen::MatrixXcd const &matrix : rep) {
    characters(element_index) = matrix.trace();
    ++element_index;
  }
  return characters;
}

/// Calculate character for all matrices in rep
Eigen::VectorXd make_characters(std::vector<Eigen::MatrixXd> const &rep) {
  Eigen::VectorXd characters(rep.size());

  Index element_index = 0;
  for (Eigen::MatrixXd const &matrix : rep) {
    characters(element_index) = matrix.trace();
    ++element_index;
  }
  return characters;
}

/// Check if approximately zero outside block along diagonal
///
/// Only checks columns and rows in range [begin, end)
bool make_is_block_diagonal(std::vector<Eigen::MatrixXcd> const &rep,
                            Index begin, Index end, double tol) {
  Index element_index = 0;
  for (Eigen::MatrixXcd const &matrix : rep) {
    // left
    if (begin != 0) {
      if (!almost_zero(matrix.block(begin, 0, end - begin, begin), tol)) {
        return false;
      }
    }
    // right
    if (end != matrix.cols()) {
      if (!almost_zero(
               matrix.block(begin, end, end - begin, matrix.cols() - end),
               tol)) {
        return false;
      }
    }
    // top
    if (begin != 0) {
      if (!almost_zero(matrix.block(0, begin, begin, end - begin), tol)) {
        return false;
      }
    }
    // bottom
    if (end != matrix.rows()) {
      if (!almost_zero(
               matrix.block(end, begin, matrix.rows() - end, end - begin),
               tol)) {
        return false;
      }
    }

    ++element_index;
  }
  return true;
}

/// Find characters for block in range [begin, end)
Eigen::VectorXcd make_irrep_characters(std::vector<Eigen::MatrixXcd> const &rep,
                                       Index begin, Index end) {
  Eigen::VectorXcd characters(rep.size());

  Index element_index = 0;
  for (Eigen::MatrixXcd const &matrix : rep) {
    characters(element_index) =
        matrix.block(begin, begin, end - begin, end - begin).trace();
    ++element_index;
  }
  return characters;
}

double make_squared_norm(Eigen::VectorXcd const &characters) {
  double squared_norm = 0.0;
  for (Index i = 0; i < characters.size(); ++i) {
    squared_norm += std::norm(characters(i));
  }
  return squared_norm;
}

double make_squared_norm(Eigen::VectorXd const &characters) {
  double squared_norm = 0.0;
  for (Index i = 0; i < characters.size(); ++i) {
    squared_norm += characters(i) * characters(i);
  }
  return squared_norm;
}

std::complex<double> frobenius_product(Eigen::MatrixXcd const &matrix) {
  return (matrix.array().conjugate() * matrix.array()).sum();
}

Eigen::MatrixXcd normalize_commuter(Eigen::MatrixXcd const &commuter) {
  return commuter / sqrt(frobenius_product(commuter).real());
}

/// Return true if space_A is extended by space_B
bool is_extended_by(Eigen::MatrixXcd const &space_A,
                    Eigen::MatrixXcd const &space_B) {
  return almost_zero((space_B.adjoint() * space_A).norm(), TOL);
}

/// Return matrix combining columns of space_A and space_B
Eigen::MatrixXcd extend(Eigen::MatrixXcd const &space_A,
                        Eigen::MatrixXcd const &space_B) {
  Eigen::MatrixXcd result(space_A.rows(), space_A.cols() + space_B.cols());
  result.leftCols(space_A.cols()) = space_A;
  result.rightCols(space_B.cols()) = space_B;
  return result;
}

/// Return matrix combining columns of space_A and space_B
Eigen::MatrixXd extend(Eigen::MatrixXd const &space_A,
                       Eigen::MatrixXd const &space_B) {
  Eigen::MatrixXd result(space_A.rows(), space_A.cols() + space_B.cols());
  result.leftCols(space_A.cols()) = space_A;
  result.rightCols(space_B.cols()) = space_B;
  return result;
}

Index get_total_dim(std::set<PossibleIrrep> const &irreps) {
  Index total_dim = 0;
  for (PossibleIrrep const &irrep : irreps) {
    total_dim += irrep.irrep_dim;
  }
  return total_dim;
}

/// \param begin: col index in `eigenvalues` and `eigenvectors` corresponding
///     to the beginning of this possible irreducible space
/// \param begin, end: col indices in `eigenvalues` and `eigenvectors`
///     corresponding to a range of equal eigenvalues
/// \param eigenvalues Eigenvalues of (K.adjoint() * M_new * K)
/// \param KV_matrix K * V, where V is the eigenvector matrix of
///     (K.adjoint() * M_new * K)
PossibleIrrep::PossibleIrrep(
    Eigen::VectorXd const &eigenvalues, Eigen::MatrixXcd const &KV_matrix,
    std::vector<Eigen::MatrixXcd> const &transformed_rep,
    Index _head_group_size, double _tol, bool allow_complex, Index _begin,
    Index _end)
    : head_group_size(_head_group_size),
      tol(_tol),
      begin(_begin),
      end(_end),
      irrep_dim(end - begin) {
  is_block_diagonal = make_is_block_diagonal(transformed_rep, begin, end, tol);
  characters = make_irrep_characters(transformed_rep, begin, end);
  characters_squared_norm = make_squared_norm(characters);
  is_irrep = is_block_diagonal && almost_equal(characters_squared_norm,
                                               double(head_group_size), tol);
  subspace = make_irrep_subspace(KV_matrix, begin, end, allow_complex);
}

/// Check if Irrep is identity
///
/// - First character is 1+0i, sum of characters == characters.size()
bool PossibleIrrep::is_identity() const {
  std::complex<double> complex_one{1., 0.};
  std::complex<double> first = characters(0);
  std::complex<double> complex_size{double(characters.size()), 0.};
  std::complex<double> sum = characters.sum();

  return almost_equal(first, complex_one, TOL) &&
         almost_equal(sum, complex_size, TOL);
}

/// Check if Irrep is gerade
///
/// - "gerade": inversion results in no sign change
/// - "ungerade": inversion results in sign change
bool PossibleIrrep::is_gerade() const {
  std::complex<double> first = characters(0);
  std::complex<double> last = characters(characters.size() - 1);
  return almost_equal(first, last, TOL);
}

bool PossibleIrrep::operator<(PossibleIrrep const &other) const {
  // Identity comes first
  bool this_is_identity = this->is_identity();
  bool other_is_identity = other.is_identity();
  if (this_is_identity != other_is_identity) {
    return this_is_identity;
  }

  // Low-dimensional irreps come before higher dimensional
  if (!almost_equal(this->characters(0), other.characters(0))) {
    return this->characters(0).real() < other.characters(0).real();
  }

  // 'gerade' irreps come before 'ungerade' irreps
  // This check may need to be improved to know whether inversion is actually
  // present
  bool this_is_gerade = this->is_gerade();
  bool other_is_gerade = other.is_gerade();
  if (this_is_gerade != other_is_gerade) {
    return this_is_gerade;
  }

  // Finally, compare lexicographically (real first, then imag)
  for (Index i = 0; i < this->characters.size(); ++i) {
    if (!almost_equal(this->characters(i).real(), other.characters(i).real()))
      return this->characters(i).real() > other.characters(i).real();
  }
  for (Index i = 0; i < this->characters.size(); ++i) {
    if (!almost_equal(this->characters(i).imag(), other.characters(i).imag()))
      return this->characters(i).imag() > other.characters(i).imag();
  }

  // Now, any possible irrep that are still equal, we break the tie by
  // comparing subspace vectors
  if (this->subspace.cols() != other.subspace.cols()) {
    return this->subspace.cols() < other.subspace.cols();
  }
  for (Index col = 0; col < this->subspace.cols(); ++col) {
    for (Index i = 0; i < this->subspace.size(); ++i) {
      if (!almost_equal(this->subspace(i, col).real(),
                        other.subspace(i, col).real()))
        return this->subspace(i, col).real() > other.subspace(i, col).real();
    }
    for (Index i = 0; i < this->subspace.size(); ++i) {
      if (!almost_equal(this->subspace(i, col).imag(),
                        other.subspace(i, col).imag()))
        return this->subspace(i, col).imag() > other.subspace(i, col).imag();
    }
  }

  throw std::runtime_error("Error comparing PossibleIrrep, tied");

  return false;
}

/// Given kernel, K, and commuter matrix, M, perform eigenvalue decomposition
///     K.adjoint() * M * K = V * D * V.inverse()
/// and construct matrix representation that acts on vectors in the K*V basis,
/// which will be block diagonalized and sorted by eigenvalue. Each block
/// corresponds to a possible irrep, which can be checked by characters value.
std::vector<PossibleIrrep> make_possible_irreps(
    Eigen::MatrixXcd const &commuter, Eigen::MatrixXcd const &kernel,
    MatrixRep const &rep, GroupIndices const &head_group, double is_irrep_tol,
    bool allow_complex) {
  // magnify the range of eigenvalues to be (I think) independent of
  // matrix dimension by multiplying by dim^{3/2}
  //
  // solve for eigenvalues and eigenvectors of:
  //    dim^(3/2) * kernel.adjoint() * M * kernel
  //
  double dim = kernel.rows();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> esolve;
  double scale = dim * sqrt(dim);
  esolve.compute(scale * kernel.adjoint() * commuter * kernel);
  Eigen::MatrixXd eigenvalues = esolve.eigenvalues();
  Eigen::MatrixXcd KV_matrix = kernel * esolve.eigenvectors();

  // Columns of KV_matrix are orthonormal eigenvectors of commuter in terms of
  // natural basis (they were calculated in terms of kernel as basis)

  // When the matrix representation is transformed to operate on coordinates
  // with KV_matrix basis, it becomes block diagonalized
  std::vector<Eigen::MatrixXcd> transformed_rep;
  transformed_rep.reserve(head_group.size());
  for (auto const &element_index : head_group) {
    transformed_rep.push_back(KV_matrix.adjoint() * rep[element_index] *
                              KV_matrix);
  }

  // make possible irreps:
  // - The possible irrep corresponds to a range eigenvectors with equal
  //   eigenvalues, could be irrep or could be reducible with degenerate
  //   eigenvalues
  // - When the possible irrep for a range of equal eigenvectors is constructed,
  //   its characters vector is constructed, and if the squared norm of the
  //   characters vectors equals the head group size, then the corresponding
  //   columns of the KV_matrix are an irrep subspace
  std::vector<PossibleIrrep> possible_irreps;
  Index begin = 0;
  do {
    Index end = find_end_of_equal_eigenvalues(begin, eigenvalues);
    possible_irreps.emplace_back(eigenvalues, KV_matrix, transformed_rep,
                                 head_group.size(), is_irrep_tol, allow_complex,
                                 begin, end);
    begin = end;
  } while (begin != eigenvalues.size());

  return possible_irreps;
}

/// Make a vector of IrrepInfo from PossibleIrreps
std::vector<IrrepInfo> make_irrep_info(std::set<PossibleIrrep> const &irreps) {
  std::vector<IrrepInfo> irrep_info;
  for (PossibleIrrep const &irrep : irreps) {
    irrep_info.emplace_back(irrep.subspace.adjoint(), irrep.characters);

    if (irrep.subspace.cols() == 2 * (irrep.end - irrep.begin))
      irrep_info.back().pseudo_irrep = true;
    else
      irrep_info.back().pseudo_irrep = false;
  }

  // set sequential indices to differentiate irreps
  // with identical character vectors
  Index irrep_index = 0;
  for (Index i = 0; i < irrep_info.size() - 1; ++i) {
    irrep_info[i].index = irrep_index;
    if (almost_equal(irrep_info[i + 1].characters, irrep_info[i].characters,
                     TOL)) {
      irrep_index++;
    } else {
      irrep_index = 0;
    }
  }

  return irrep_info;
}

/// \brief Transforms IrrepInfo constructed for a subspace to be IrrepInfo
/// appropriate for the full space (full space dimension == subspace.rows())
///
/// \param irrep IrrepInfo constructed for a subspace (irrep.trans_mat shape is
///     (subspace.cols() x subspace.rows())
/// \param subspace The subspace that subspace_irrep was constructed for
///
/// \result IrrepInfo constructed for the full space (result.trans_mat shape is
///     (subspace.cols() x subspace.cols())
///
IrrepInfo subspace_to_full_space(IrrepInfo const &subspace_irrep,
                                 Eigen::MatrixXd const &subspace) {
  IrrepInfo result(subspace_irrep);

  result.trans_mat = subspace_irrep.trans_mat *
                     subspace.adjoint().template cast<std::complex<double>>();

  result.directions.clear();
  for (const auto &direction_orbit : subspace_irrep.directions) {
    std::vector<Eigen::VectorXd> new_orbit;
    new_orbit.reserve(direction_orbit.size());
    for (const auto &directions : direction_orbit) {
      new_orbit.push_back(subspace * directions);
    }
    result.directions.push_back(std::move(new_orbit));
  }
  return result;
}

/// Check if a representation is irreducible
///
/// A representation is irreducible if the squared norm of the characters
/// equals the group size
bool is_irrep(MatrixRep const &rep, GroupIndices const &head_group) {
  double characters_squared_norm = 0;
  for (Index element_index : head_group) {
    double character = rep[element_index].trace();
    characters_squared_norm += character * character;
  }
  return almost_equal(characters_squared_norm, double(head_group.size()), TOL);
}

/// IrrepDecomposition proceeds by constructing "commuters", M_k, which commute
/// (M_k * R(r) = R(r) * M_k) with all of the matrix representations, R(r), of
/// the group. The commuters are constructed to reveal irreducible vector
/// spaces (via application of a Reynolds operator), and be orthonormal to
/// existing commuters (via Gram-Shmidt). The commuters are found via a
/// process which constructs a candidate commuter which is either the Zero
/// matrix, and then skipped, or else it is a useful non-zero commuter which
/// will block diagonalize :
///
///     M_candidate(i,j,phase) = sum_r R(r) * M_init * R(r).transpose()
///     M_init = phase * K.col(i) * K.col(j).adjoint() +
///              std::conj(phase) * K.col(j) * K.col(i).adjoint()
///
/// where:
/// - K: kernel matrix, the null space of the already found irreducible vector
/// spaces. The kernel matrix is initialized as a full rank matrix and over the
/// course of the IrrepDecomposition the kernel shrinks as the irreducible
/// vector spaces are found.
/// - candidate commuting matrices, M_candidate, are built from the outer
/// product of two columns, i, and j, of the kernel matrix, and a complex phase
/// parameter (1 or i)
/// - R(r): is the matrix representation for element r of the head_group
/// - M_k: previously found commuting matrices
///
/// Once a new non-zero commuter is found, possible irreducible subspaces are
/// found and checked. A possible irreducible subspace is each
/// K*V_equal_eigenvalue_set[i], where V_equal_eigenvalue_set[i] is the vector
/// space corresponding to eigenvectors of (K.adjoint() * M_new * K) with equal
/// eigenvalues.
///
/// Eigenvalue decomposition:
///
///     K.adjoint() * M_new * K = V * D * V.inverse()
///
/// - D: diagonal matrix of sorted eigenvalues (size K.cols() x K.cols())
/// - V: eigenvector column matrix (size K.cols() x K.cols())
///
/// For each set of equal eigenvalues, a PossibleIrrep is constructed that
/// stores:
/// - `begin`: the column of the first of the set of equal eigenvalues
/// - `irrep_dim`: the subspace dimension / the number of equal eigenvalues
/// - `subspace`: irreducible subspace, (dim x irrep_dim (?) matrix):
///
///   First step, find vector space:
///   - If allow_complex: subspace = X
///   - If !allow_complex: subspace =
///         [sqrt(2.0) * X.real(), sqrt(2.0) *  X.imag()],
///     where X = (K * V).block(0, begin,
///                             K.rows(), irrep_dim)
///   Second step, orthogonalize via QR decomposition
///
/// - `characters`: Vector of complex characters. The character of a matrix
///   representation is the trace of the representation.
/// - `characters_squared_norm`: For an irreducible representation, the squared
///   norm of the characters vector is equal the size of the group.
/// - `symmetrizer`: A pair with, symmetrizer.first being a MatrixXcd, which
///   defines a rotation of the irreducible subspace that aligns its components
///   along high-symmetry directions, and symmetrizer.second being a vector of
///   orbits of high-symmetry directions
///
/// For each PossibleIrrep, check if it extends adapted_subspace. If it
/// does, then symmetrize and save the irrep. For all new irreps, extend
/// adapted_subspace to include the irrep's subspace. Then once all the new
/// irreps are added, recalculate the kernel matrix and begin the loop again,
/// until the adapted_subspace is full rank.

/// Finds irreducible subspaces that comprise an underlying subspace
///
/// This method does not rely on the character table, but instead utilizes a
/// brute-force approach. It is not guaranteed to find all irreps, so the
/// resulting irreps should be checked if they span the entire space
/// represented by `rep`. This can be done by checking if
///     `full_trans_mat(result).adjoint().rows() == rep[i].rows()`).
/// This method does not align the irrep subspace axes along high symmetry
/// directions.
///
/// \param rep Matrix representation of head_group, this defines group action
/// on the underlying vector space
/// \param head_group Group for which the irreps are to be found
/// \param allow_complex If true, irreducible space basis vectors may be
///     complex-valued. If false, complex irreps are combined to form real
///     representations
///
/// \result vector of IrrepInfo objects. Irreps are ordered by dimension, with
///     identity first (if present).  Repeated irreps (with equal character
///     vectors) are sequential, and are distinguished by IrrepInfo::index.
///
std::vector<IrrepInfo> irrep_decomposition(MatrixRep const &rep,
                                           GroupIndices const &head_group,
                                           bool allow_complex) {
  if (!rep.size()) {
    return std::vector<IrrepInfo>();
  }
  int dim = rep[0].rows();

  // This method iteratively finds irreducible spaces, which are used to extend
  // the "adapted_subspace" (combined space of found irreducible spaces). The
  // "adapted_subspace" is not aligned along high symmetry directions by this
  // function. When the "adapted_subspace" is of dimenions equal to `dim`, all
  // irreps have been found.

  // start with all kernel, end with all adapted_subspace
  Eigen::MatrixXcd kernel = complex_I(dim, dim);
  Eigen::MatrixXcd adapted_subspace{dim, 0};

  // In this set, as they are discovered we will save PossibleIrrep that:
  // - i) actually are irreducible,
  // - and ii) have distinct subspaces (detected when their subspace extends
  //   the adapted_subspace space) BP: not necessary?
  std::set<PossibleIrrep> irreps;

  // count over possible commuter matrices for this kernel:
  // - kernel column pairs (i,j), j>=i & phase = [1, i]; skips i==j if phase==i
  CommuterParamsCounter commuter_params;
  commuter_params.reset(kernel);

  double is_irrep_tol = TOL;

  do {  // while adapated_subspace.cols() != dim

    if (!commuter_params.valid()) {
      // The commuter construction method does not currently guarantee that
      // all irreps will be revealed. The caller may have a way to handle this
      // and so this does not throw an exception.
      break;
    }

    // make next commuter, M, and check if not zero
    Eigen::MatrixXcd commuter =
        make_commuter(commuter_params, rep, head_group, kernel);

    if (almost_equal(frobenius_product(commuter).real(), 0., TOL)) {
      commuter_params.increment();
      continue;
    }

    // make possible irreps:
    //
    // Given kernel, K, and commuter matrix, M, perform eigenvalue decomposition
    //     K.adjoint() * M * K = V * D * V.inverse()
    // and construct matrix representation that acts on vectors in the K*V
    // basis, which will be block diagonalized and sorted by eigenvalue. Each
    // block corresponds to a possible irrep, which can be checked by its
    // characters. The columns in K*V corresponding to an irrep are the irrep
    // subspace.
    std::vector<PossibleIrrep> possible_irreps = make_possible_irreps(
        commuter, kernel, rep, head_group, is_irrep_tol, allow_complex);

    // save any possible irrep that:
    // - i) is an irrep,
    // - and ii) extends the adapted_subspace space (BP: not necessary?)
    bool any_new_irreps = false;
    for (auto const &possible_irrep : possible_irreps) {
      if (possible_irrep.is_irrep &&
          is_extended_by(adapted_subspace, possible_irrep.subspace)) {
        irreps.insert(possible_irrep);
        adapted_subspace = extend(adapted_subspace, possible_irrep.subspace);
        any_new_irreps = true;
      }
    }

    // if any new irreps were found, recalculate kernel, reset commuter params
    // counter, and go again
    if (any_new_irreps && adapted_subspace.cols() != dim) {
      kernel = make_kernel(adapted_subspace);
      commuter_params.reset(kernel);
      if (kernel.cols() + adapted_subspace.cols() != adapted_subspace.rows()) {
        throw std::runtime_error(
            "Unknown error finding irreps: dimension mismatch");
      }
    } else {
      commuter_params.increment();
    }
  } while (adapted_subspace.cols() != dim);

  // Make irrep info (no directions yet, orthogonalized but not aligned along
  // high symmetry directions)
  std::vector<IrrepInfo> irrep_info = make_irrep_info(irreps);

  return irrep_info;
}

/// Convert irreps generated for a subspace to full space dimension
///
/// \param subspace_irreps Irreducible spaces in the subspace
/// (subspace_irreps[i].trans_mat.rows() == subspace dimension,
/// subspace_irreps[i].trans_mat.cols() == fullspace dimensino) \param subspace
/// Basis for a subspace (subspace.rows() == fullspace
///     dimension, subspace.cols() == subspace dimension)
std::vector<IrrepInfo> make_fullspace_irreps(
    std::vector<IrrepInfo> const &subspace_irreps,
    Eigen::MatrixXd const &subspace) {
  std::vector<IrrepInfo> fullspace_irreps;
  fullspace_irreps.reserve(subspace_irreps.size());
  for (auto const &irrep : subspace_irreps) {
    fullspace_irreps.push_back(subspace_to_full_space(irrep, subspace));
  }
  return fullspace_irreps;
}

/// Expand subspace by application of group, and orthogonalize
Eigen::MatrixXd make_invariant_space(MatrixRep const &rep,
                                     GroupIndices const &head_group,
                                     Eigen::MatrixXd const &subspace) {
  if (!subspace.isIdentity()) {
    Eigen::MatrixXd symspace(subspace.rows(),
                             subspace.cols() * head_group.size());
    Index l = 0;
    for (Index element_index : head_group) {
      symspace.block(0, l, subspace.rows(), subspace.cols()) =
          rep[element_index] * subspace;
      l += subspace.cols();
    }
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> colqr(symspace);
    colqr.setThreshold(TOL);
    Eigen::MatrixXd Q = colqr.householderQ();
    Eigen::MatrixXd result = Q.leftCols(colqr.rank());
    return result;
  }
  return subspace;
}

// Create `subspace_rep`, a transformed copy of `fullspace_rep` that acts
// on coordinates with `subspace` columns as a basis. Matrices in
// `subspace_rep` are shape (subspace.cols() x subspace.cols())
MatrixRep make_subspace_rep(MatrixRep const &fullspace_rep,
                            Eigen::MatrixXd const &subspace) {
  Eigen::MatrixXd trans_mat = subspace.transpose();
  Eigen::MatrixXd rightmat =
      subspace.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
          .solve(Eigen::MatrixXd::Identity(trans_mat.cols(), trans_mat.cols()))
          .transpose();
  MatrixRep subspace_rep;
  for (Index i = 0; i < fullspace_rep.size(); ++i) {
    subspace_rep.push_back(trans_mat * fullspace_rep[i] * rightmat);
  }
  return subspace_rep;
}

/// \brief Symmetrize IrrepInfo, by finding high symmetry directions and
/// aligning the irrep subspace basis with those directions
std::vector<IrrepInfo> symmetrize_irreps(
    MatrixRep const &subspace_rep, GroupIndices const &head_group,
    std::vector<IrrepInfo> const &irreps,
    GroupIndicesOrbitVector const &cyclic_subgroups,
    GroupIndicesOrbitVector const &all_subgroups) {
  std::vector<IrrepInfo> symmetrized_irreps;
  double vec_compare_tol = TOL;
  bool use_all_subgroups = false;
  for (const auto &irrep : irreps) {
    Eigen::MatrixXcd irrep_subspace = irrep.trans_mat.adjoint();

    multivector<Eigen::VectorXcd>::X<2> irrep_special_directions =
        make_irrep_special_directions(subspace_rep, head_group, irrep_subspace,
                                      vec_compare_tol, cyclic_subgroups,
                                      all_subgroups, use_all_subgroups);

    Eigen::MatrixXcd symmetrizer_matrix = make_irrep_symmetrizer_matrix(
        irrep_special_directions, irrep_subspace, vec_compare_tol);

    IrrepInfo symmetrized_irrep{irrep};
    symmetrized_irrep.trans_mat =
        (irrep_subspace * symmetrizer_matrix).adjoint();
    symmetrized_irrep.directions = to_real(irrep_special_directions);
    symmetrized_irreps.push_back(symmetrized_irrep);
  }
  return symmetrized_irreps;
}

}  // namespace IrrepDecompositionImpl

}  // namespace SymRepTools_v2

}  // namespace CASM
