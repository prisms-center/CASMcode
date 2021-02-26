#include "casm/symmetry/Symmetrizer.hh"

#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/symmetry/SimpleOrbit_impl.hh"

namespace CASM {

namespace SymRepTools_v2 {

/// Used for generating SimpleOrbit of high symmetry direction vectors
class VectorInvariants {
 public:
  VectorInvariants(Eigen::VectorXcd const &vector);

  double cols() const;
  double norm() const;

 private:
  double m_cols;
  double m_norm;
};

VectorInvariants::VectorInvariants(Eigen::VectorXcd const &vector)
    : m_cols(vector.cols()), m_norm(vector.norm()) {}

double VectorInvariants::cols() const { return m_cols; }

double VectorInvariants::norm() const { return m_norm; }

}  // namespace SymRepTools_v2

bool almost_equal(SymRepTools_v2::VectorInvariants const &A_invariants,
                  SymRepTools_v2::VectorInvariants const &B_invariants,
                  double tol);

bool compare(SymRepTools_v2::VectorInvariants const &A_invariants,
             SymRepTools_v2::VectorInvariants const &B_invariants, double tol);

bool almost_equal(SymRepTools_v2::VectorInvariants const &A_invariants,
                  SymRepTools_v2::VectorInvariants const &B_invariants,
                  double tol) {
  return A_invariants.cols() == B_invariants.cols() &&
         almost_equal(A_invariants.norm(), B_invariants.norm(), tol);
}

bool compare(SymRepTools_v2::VectorInvariants const &A_invariants,
             SymRepTools_v2::VectorInvariants const &B_invariants, double tol) {
  if (A_invariants.cols() == B_invariants.cols()) {
    return CASM::compare(A_invariants.norm(), B_invariants.norm(), tol);
  }
  return A_invariants.cols() < B_invariants.cols();
}

namespace SymRepTools_v2 {

/// Used for constructing SimpleOrbit of high symmetry direction vectors
struct VectorSymCompare {
  VectorSymCompare(MatrixRep const &matrix_rep, double tol);

  typedef Index SymOpRepType;
  typedef Eigen::VectorXcd Element;
  typedef VectorInvariants InvariantsType;

  // lexicographical comparison (reversed, uses vector_B < vector_A)
  bool compare(Eigen::VectorXcd const &vector_A,
               Eigen::VectorXcd const &vector_B) const;

  // return VectorInvariants
  VectorInvariants make_invariants(Eigen::VectorXcd const &vector) const;

  // apply matrix rep to vector
  Eigen::VectorXcd copy_apply(Index const &op_index,
                              Eigen::VectorXcd vector) const;

  // no change needed to prepare for comparison
  Eigen::VectorXcd prepare(Eigen::VectorXcd vector) const;

  // compare orbits by first comparing invariants, then orbit prototypes
  bool inter_orbit_compare(Eigen::VectorXcd const &A_prototype,
                           VectorInvariants const &A_invariants,
                           Eigen::VectorXcd const &B_prototype,
                           VectorInvariants const &B_invariants) const;

 private:
  MatrixRep const &m_matrix_rep;
  double m_tol;
};

VectorSymCompare::VectorSymCompare(MatrixRep const &matrix_rep, double tol)
    : m_matrix_rep(matrix_rep), m_tol(tol) {}

// lexicographical comparison (reversed, uses vector_B < vector_A)
bool VectorSymCompare::compare(Eigen::VectorXcd const &vector_A,
                               Eigen::VectorXcd const &vector_B) const {
  return colmajor_lex_compare(vector_B, vector_A, m_tol);
}

// return VectorInvariants
VectorInvariants VectorSymCompare::make_invariants(
    Eigen::VectorXcd const &vector) const {
  return VectorInvariants{vector};
}

// apply matrix rep to vector
Eigen::VectorXcd VectorSymCompare::copy_apply(Index const &op_index,
                                              Eigen::VectorXcd vector) const {
  return m_matrix_rep[op_index] * vector;
}

// no change needed to prepare for comparison
Eigen::VectorXcd VectorSymCompare::prepare(Eigen::VectorXcd vector) const {
  return vector;
}

// compare orbits by first comparing invariants, then orbit prototypes
bool VectorSymCompare::inter_orbit_compare(
    Eigen::VectorXcd const &A_prototype, VectorInvariants const &A_invariants,
    Eigen::VectorXcd const &B_prototype,
    VectorInvariants const &B_invariants) const {
  // first compare invariants
  if (CASM::compare(A_invariants, B_invariants, m_tol)) {
    return true;
  }
  if (CASM::compare(B_invariants, A_invariants, m_tol)) {
    return false;
  }

  // next compare A and B
  return this->compare(A_prototype, B_prototype);
}

/// Vector space preparation
///
/// - Attempts to find a sparse set of spanning vectors, and sort them so that
/// subspace matrix is nearly upper triangular (if possible)
/// - Also enures that first nonzero element of each row (if there is one) is
/// positive
Eigen::MatrixXcd vector_space_prepare(Eigen::MatrixXcd const &vector_space,
                                      double _tol) {
  Eigen::MatrixXcd reduced_column_echelon_form =
      CASM::reduced_column_echelon(vector_space, _tol);
  Eigen::MatrixXcd Q =
      reduced_column_echelon_form.householderQr().householderQ();
  Eigen::MatrixXcd result = Q.leftCols(vector_space.cols());
  CASM::Index col = 0;
  for (CASM::Index row = 0; row < result.rows(); ++row) {
    CASM::Index i = 0;
    for (i = col; i < result.cols(); ++i) {
      if (!CASM::almost_zero(result(row, i), _tol)) {
        result.col(i) *= std::abs(result(row, i)) / result(row, i);
        ++col;
        break;
      }
    }
  }
  return result;
}

/// Find high-symmetry directions in a irreducible space
///
/// \param rep Matrix representation of head_group, this defines group action
/// on the underlying vector space
/// \param head_group Group for which the irreps are to be found
/// \param subspace A column vector matrix representing a basis of the
///     irreducible space in which high symmetry directions will be found. The
///     number of rows must equal `rep.dim()`, the number of columns is equal to
///     the dimension of the irreducible space.
/// \param vec_compare_tol Tolerance for elementwise floating-point comparisons
///     of vectors
/// \param all_subgroups Denotes whether all subgroups of head_group should be
///     used for symmetry analysis (if true), or only cyclic subgroups (if
///     false). Cyclic subgroups are those found by taking a group element and
///     multiplying it by itself until a group is generated.
///
/// \result Set of directions in the vector space on which 'rep' is defined,
/// such that each direction is invariant to a unique subgroup of 'head_group'
/// (i.e., no other direction in the space, except the negative of that
/// direction, is invariant to that subgroup). The value `result[i]` is an
/// orbit of symmetrically equivalent directions, and the value `result[i][j]`
/// is an individual direction (Eigen::VectorXcd). Direction vectors are
/// normalized to unit length. The total set of all directions is guaranteed to
/// span the space.
///
/// \throws if `rep` is not an irreducible representation
///
multivector<Eigen::VectorXcd>::X<2> make_irrep_special_directions(
    MatrixRep const &rep, GroupIndices const &head_group,
    Eigen::MatrixXcd const &irrep_subspace, double vec_compare_tol,
    GroupIndicesOrbitVector const &cyclic_subgroups,
    GroupIndicesOrbitVector const &all_subgroups, bool use_all_subgroups) {
  auto const &sgroups = (use_all_subgroups ? all_subgroups : cyclic_subgroups);

  std::vector<Eigen::VectorXcd> tdirs;
  Eigen::MatrixXd R;
  Index dim = rep[0].rows();

  // Loop over small (i.e., cyclic) subgroups and hope that each special
  // direction is invariant to at least one small subgroup
  for (auto const &orbit : sgroups) {
    // Reynolds for small subgroup *(orbit.begin()) in irrep_subspace i
    R.setZero(dim, dim);

    for (Index element_index : *(orbit.begin())) {
      R += rep[element_index];
    }

    if ((R * irrep_subspace).norm() < TOL) continue;

    // Find spanning vectors of column space of R*irrep_space, which is
    // projection of irrep_space into its invariant component
    auto QR = (R * irrep_subspace).colPivHouseholderQr();
    QR.setThreshold(TOL);
    // If only one spanning vector, it is special direction
    if (QR.rank() > 1) continue;
    Eigen::MatrixXcd Q = QR.matrixQ();

    // Convert from irrep_subspace back to total space and push_back
    tdirs.push_back(Q.col(0));
    tdirs.push_back(-Q.col(0));
  }

  // t_result may contain duplicates, or elements that are equivalent by
  // symmetry. To discern more info, we need to exclude duplicates and find
  // the orbit of the directions. this should also
  // reveal the invariant subgroups.

  VectorSymCompare sym_compare{rep, vec_compare_tol};
  std::set<SimpleOrbit<VectorSymCompare>> orbit_result;
  for (Eigen::VectorXcd const &direction : tdirs) {
    orbit_result.emplace(direction, head_group.begin(), head_group.end(),
                         sym_compare);
  }
  multivector<Eigen::VectorXcd>::X<2> result;
  for (auto const &orbit : orbit_result) {
    result.emplace_back(orbit.begin(), orbit.end());
  }

  if (use_all_subgroups || result.size() >= irrep_subspace.cols()) {
    return result;
  } else {
    return make_irrep_special_directions(rep, head_group, irrep_subspace,
                                         vec_compare_tol, cyclic_subgroups,
                                         all_subgroups, true);
  }
}

/// Make an irreducible space symmetrizer matrix using special directions
Eigen::MatrixXcd make_irrep_symmetrizer_matrix(
    multivector<Eigen::VectorXcd>::X<2> const &irrep_special_directions,
    Eigen::MatrixXcd const &irrep_subspace, double vec_compare_tol) {
  // Four strategies, in order of desparation
  // 1) find a spanning set of orthogonal axes within a single orbit of special
  // directions 2) find a spanning set of orthogonal axes within the total set
  // of special directions 3) perform qr decomposition on lowest-multiplicity
  // orbit to find spanning set of axes 4) find column-echelon form of
  // irrep_subspace matrix to get a sparse/pretty set of axes
  Eigen::MatrixXcd result;
  Index dim = irrep_subspace.cols();
  Index min_mult = 10000;
  bool orb_orthog = false;
  bool tot_orthog = false;

  Eigen::MatrixXcd axes, orb_axes, tot_axes;
  Index tot_col(0);
  tot_axes.setZero(irrep_subspace.rows(), dim);

  for (auto const &orbit : irrep_special_directions) {
    // Strategy 1
    if ((orb_orthog && orbit.size() < min_mult) || !orb_orthog) {
      orb_axes.setZero(irrep_subspace.rows(), dim);
      Index col = 0;
      for (auto const &el : orbit) {
        if (almost_zero((el.adjoint() * orb_axes).eval(), vec_compare_tol)) {
          if (col < orb_axes.cols())
            orb_axes.col(col++) = el;
          else {
            std::stringstream errstr;
            errstr
                << "Error in irrep_symmtrizer_from_directions(). Constructing "
                   "coordinate axes from special directions of space spanned "
                   "by row vectors:\n "
                << irrep_subspace.transpose()
                << "\nAxes collected thus far are the row vectors:\n"
                << orb_axes.transpose()
                << "\nBut an additional orthogonal row vector has been found:\n"
                << el.transpose()
                << "\nindicating that irrep_subspace matrix is malformed.";
            throw std::runtime_error(errstr.str());
          }
        }
      }
      if (col == dim) {
        // std::cout << "PLAN A-- col: " << col << "; dim: " << dim << ";
        // min_mult: " << min_mult << "; orthog: " << orb_orthog << ";\naxes:
        // \n"
        // << axes << "\norb_axes: \n" << orb_axes << "\n\n";
        orb_orthog = true;
        min_mult = orbit.size();
        axes = orb_axes;
      }
    }

    // Greedy(ish) implementation of strategy 2 -- may not find a solution, even
    // if it exists
    if (!orb_orthog && !tot_orthog) {
      for (auto const &el : orbit) {
        if (almost_zero((el.adjoint() * tot_axes).eval(), vec_compare_tol)) {
          if (tot_col < tot_axes.cols())
            tot_axes.col(tot_col++) = el;
          else {
            std::stringstream errstr;
            errstr
                << "Error in irrep_symmtrizer_from_directions(). Constructing "
                   "coordinate axes from special directions of space spanned "
                   "by row vectors:\n "
                << irrep_subspace.transpose()
                << "\nAxes collected thus far are the row vectors:\n"
                << tot_axes.transpose()
                << "\nBut an additional orthogonal row vector has been found:\n"
                << el.transpose()
                << "\nindicating that irrep_subspace matrix is malformed.";
            throw std::runtime_error(errstr.str());
          }
        }
      }
      if (tot_col == dim) {
        // std::cout << "PLAN B-- col: " << tot_col << "; dim: " << dim << ";
        // min_mult: " << min_mult << "; orthog: " << tot_orthog << ";\naxes:
        // \n"
        // << axes << "\ntot_axes: \n" << tot_axes << "\n\n";
        tot_orthog = true;
        axes = tot_axes;
      }
    }

    // Strategy 3
    if (!orb_orthog && !tot_orthog && orbit.size() < min_mult) {
      orb_axes.setZero(irrep_subspace.rows(), orbit.size());
      for (Index col = 0; col < orbit.size(); ++col) {
        orb_axes.col(col) = orbit[col];
      }
      // std::cout << "PLAN C--  dim: " << dim << "; min_mult: " << min_mult <<
      // "; orthog: " << orb_orthog << ";\naxes: \n" << axes;
      min_mult = orbit.size();
      axes = Eigen::MatrixXcd(orb_axes.colPivHouseholderQr().matrixQ())
                 .leftCols(dim);
      // std::cout << "\norb_axes: \n" << orb_axes << "\n\n";
    }
  }
  // std::cout << "axes: \n" << axes << "\n";
  if (axes.cols() == 0) {
    result = irrep_subspace.colPivHouseholderQr().solve(
        vector_space_prepare(irrep_subspace, vec_compare_tol));
  } else {
    // Strategy 4
    result = irrep_subspace.colPivHouseholderQr().solve(axes);
  }

  // std::cout << "result: \n" << result << "\n"
  //<< "irrep_subspace*result: \n" << irrep_subspace *result << "\n";
  return result;  //.transpose();
}

}  // namespace SymRepTools_v2

}  // namespace CASM
