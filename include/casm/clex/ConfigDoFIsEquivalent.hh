#ifndef CASM_ConfigDoFIsEquivalent
#define CASM_ConfigDoFIsEquivalent

#include "casm/clexulator/ConfigDoFValues.hh"
#include "casm/clexulator/ConfigDoFValuesTools.hh"
#include "casm/symmetry/PermuteIterator.hh"

namespace CASM {

/** \ingroup ConfigIsEquivalent
 *  @{
 */

/// Namespace containing DoF comparison functors
namespace ConfigDoFIsEquivalent {

/// Compare isotropic occupation values
///
/// - The protected '_check' method provides for both checking equality and if
///   not equivalent, storing the 'less than' result
class Occupation {
 public:
  Occupation(Eigen::VectorXi const &_occupation)
      : m_occupation_ptr(&_occupation) {}

  /// \brief Return config == other, store config < other
  bool operator()(Eigen::VectorXi const &other) const {
    return _for_each([&](Index i) { return (*m_occupation_ptr)[i]; },
                     [&](Index i) { return other[i]; });
  }

  /// \brief Return config == A*config, store config < A*config
  bool operator()(PermuteIterator const &A) const {
    return _for_each(
        [&](Index i) { return (*m_occupation_ptr)[i]; },
        [&](Index i) { return (*m_occupation_ptr)[A.permute_ind(i)]; });
  }

  /// \brief Return A*config == B*config, store A*config < B*config
  bool operator()(PermuteIterator const &A, PermuteIterator const &B) const {
    return _for_each(
        [&](Index i) { return (*m_occupation_ptr)[A.permute_ind(i)]; },
        [&](Index i) { return (*m_occupation_ptr)[B.permute_ind(i)]; });
  }

  /// \brief Return config == A*other, store config < A*other
  bool operator()(PermuteIterator const &A,
                  Eigen::VectorXi const &other) const {
    return _for_each([&](Index i) { return (*m_occupation_ptr)[i]; },
                     [&](Index i) { return other[A.permute_ind(i)]; });
  }

  /// \brief Return A*config == B*other, store A*config < B*other
  bool operator()(PermuteIterator const &A, PermuteIterator const &B,
                  Eigen::VectorXi const &other) const {
    return _for_each(
        [&](Index i) { return (*m_occupation_ptr)[A.permute_ind(i)]; },
        [&](Index i) { return other[B.permute_ind(i)]; });
  }

  /// \brief Returns less than comparison
  ///
  /// - Only valid after call operator returns false
  bool is_less() const { return m_less; }

 protected:
  template <typename F, typename G>
  bool _for_each(F f, G g) const {
    Index i;
    for (i = 0; i < m_occupation_ptr->size(); i++) {
      if (!_check(f(i), g(i))) {
        return false;
      }
    }
    return true;
  }

 private:
  template <typename T>
  bool _check(const T &A, const T &B) const {
    if (A == B) {
      return true;
    }
    m_less = (A < B);
    return false;
  }

  Eigen::VectorXi const *m_occupation_ptr;

  /// Stores (A < B) if A != B
  mutable bool m_less;
};

/// Compare anisotropic occupation values
///
/// - The protected '_check' method provides for both checking equality and if
///   not equivalent, storing the 'less than' result
///
/// Method:
/// - To improve efficiency when comparisons are being made repeatedly under
///   transformation by the same factor group operation but different
///   translations, temporary vectors 'm_new_occ_A' and 'm_new_occ_B' store the
///   occupation values under transformation by factor group operation only,
///   not site permutation. The factor group operation used last is stored in
///   'm_fg_index_A' and 'm_fg_index_B'. The 'm_tmp_valid' is set to false if
///   comparison is made against an "other" ConfigDoF to force update of the
///   transformed variables in the temporary vectors the next time the functor
///   is called because it cannot be guaranteed that the "other" is the same.
class AnisoOccupation {
 public:
  AnisoOccupation(Eigen::VectorXi const &_occupation, Index n_sublat)
      : m_n_sublat(n_sublat),
        m_n_vol(_occupation.size() / m_n_sublat),
        m_occupation_ptr(&_occupation),
        m_tmp_valid(true),
        m_fg_index_A(0),
        m_new_occ_A(_occupation),
        m_fg_index_B(0),
        m_new_occ_B(_occupation) {}

  /// \brief Return config == other, store config < other
  bool operator()(Eigen::VectorXi const &other) const {
    return _for_each([&](Index i) { return (*m_occupation_ptr)[i]; },
                     [&](Index i) { return other[i]; });
  }

  /// \brief Return config == B*config, store config < B*config
  bool operator()(PermuteIterator const &B) const {
    _update_B(B, *m_occupation_ptr);
    m_tmp_valid = true;

    return _for_each(
        [&](Index i) { return (*m_occupation_ptr)[i]; },
        [&](Index i) { return this->m_new_occ_B[B.permute_ind(i)]; });
  }

  /// \brief Return A*config == B*config, store A*config < B*config
  bool operator()(PermuteIterator const &A, PermuteIterator const &B) const {
    _update_A(A, *m_occupation_ptr);
    _update_B(B, *m_occupation_ptr);
    m_tmp_valid = true;

    return _for_each(
        [&](Index i) { return this->m_new_occ_A[A.permute_ind(i)]; },
        [&](Index i) { return this->m_new_occ_B[B.permute_ind(i)]; });
  }

  /// \brief Return config == B*other, store config < B*other
  bool operator()(PermuteIterator const &B,
                  Eigen::VectorXi const &other) const {
    _update_B(B, other);
    m_tmp_valid = false;

    return _for_each(
        [&](Index i) { return (*m_occupation_ptr)[i]; },
        [&](Index i) { return this->m_new_occ_B[B.permute_ind(i)]; });
  }

  /// \brief Return A*config == B*other, store A*config < B*other
  bool operator()(PermuteIterator const &A, PermuteIterator const &B,
                  Eigen::VectorXi const &other) const {
    _update_A(A, *m_occupation_ptr);
    _update_B(B, other);
    m_tmp_valid = false;

    return _for_each(
        [&](Index i) { return this->m_new_occ_A[A.permute_ind(i)]; },
        [&](Index i) { return this->m_new_occ_B[B.permute_ind(i)]; });
  }

  /// \brief Returns less than comparison
  ///
  /// - Only valid after call operator returns false
  bool is_less() const { return m_less; }

 protected:
  template <typename F, typename G>
  bool _for_each(F f, G g) const {
    Index i;
    for (i = 0; i < m_occupation_ptr->size(); i++) {
      if (!_check(f(i), g(i))) {
        return false;
      }
    }
    return true;
  }

  void _update_A(PermuteIterator const &A,
                 Eigen::VectorXi const &before) const {
    if (A.factor_group_index() != m_fg_index_A || !m_tmp_valid) {
      m_fg_index_A = A.factor_group_index();
      Index l = 0;
      for (Index b = 0; b < m_n_sublat; ++b) {
        for (Index n = 0; n < m_n_vol; ++n, ++l) {
          m_new_occ_A[l] = (*(A.occ_rep(b).permutation()))[before[l]];
        }
      }
    }
  }

  void _update_B(PermuteIterator const &B,
                 Eigen::VectorXi const &before) const {
    if (B.factor_group_index() != m_fg_index_B || !m_tmp_valid) {
      m_fg_index_B = B.factor_group_index();
      Index l = 0;
      for (Index b = 0; b < m_n_sublat; ++b) {
        for (Index n = 0; n < m_n_vol; ++n, ++l) {
          m_new_occ_B[l] = (*(B.occ_rep(b).permutation()))[before[l]];
        }
      }
    }
  }

 private:
  template <typename T>
  bool _check(const T &A, const T &B) const {
    if (A == B) {
      return true;
    }
    m_less = (A < B);
    return false;
  }

  Index m_n_sublat;

  Index m_n_vol;

  // Points to ConfigDoF this was constructed with
  Eigen::VectorXi const *m_occupation_ptr;

  // Set to false when comparison is made to "other" ConfigDoF, to force update
  // of temporary dof during the next comparison
  mutable bool m_tmp_valid;

  // Store temporary DoFValues under fg operation only:

  mutable Index m_fg_index_A;
  mutable Eigen::VectorXi m_new_occ_A;

  mutable Index m_fg_index_B;
  mutable Eigen::VectorXi m_new_occ_B;

  /// Stores (A < B) if A != B
  mutable bool m_less;
};

/// Compare continuous site DoF values
///
/// Method:
/// - To improve efficiency when comparisons are being made repeatedly under
///   transformation by the same factor group operation but different
///   translations, temporary vectors 'm_new_dof_A' and 'm_new_dof_B' store the
///   DoF values under transformation by factor group operation only,
///   not site permutation. The factor group operation used last is stored in
///   'm_fg_index_A' and 'm_fg_index_B'. The 'm_tmp_valid' is set to false if
///   comparison is made against an "other" ConfigDoF to force update of the
///   transformed variables in the temporary vectors the next time the functor
///   is called because it cannot be guaranteed that the "other" is the same.
class Local {
 public:
  Local(Eigen::MatrixXd const &_values, DoFKey const &_key, Index n_sublat,
        double _tol)
      : m_values_ptr(&_values),
        m_key(_key),
        m_n_sublat(n_sublat),
        m_n_vol(_values.cols() / n_sublat),
        m_tol(_tol),
        m_tmp_valid(true),
        m_fg_index_A(0),
        m_new_dof_A(*m_values_ptr),
        m_fg_index_B(0),
        m_new_dof_B(*m_values_ptr) {}

  /// \brief Return config == other, store config < other
  bool operator()(Eigen::MatrixXd const &other) const {
    return _for_each([&](Index i, Index j) { return this->_values()(i, j); },
                     [&](Index i, Index j) { return other(i, j); });
  }

  /// \brief Return config == B*config, store config < B*config
  bool operator()(PermuteIterator const &B) const {
    _update_B(B, _values());
    m_tmp_valid = true;

    return _for_each(
        [&](Index i, Index j) { return this->_values()(i, j); },
        [&](Index i, Index j) { return this->new_dof_B(i, B.permute_ind(j)); });
  }

  /// \brief Return A*config == B*config, store A*config < B*config
  bool operator()(PermuteIterator const &A, PermuteIterator const &B) const {
    _update_A(A, _values());
    _update_B(B, _values());
    m_tmp_valid = true;
    return _for_each(
        [&](Index i, Index j) { return this->new_dof_A(i, A.permute_ind(j)); },
        [&](Index i, Index j) { return this->new_dof_B(i, B.permute_ind(j)); });
  }

  /// \brief Return config == B*other, store config < B*other
  bool operator()(PermuteIterator const &B,
                  Eigen::MatrixXd const &other) const {
    _update_B(B, other);
    m_tmp_valid = false;

    return _for_each(
        [&](Index i, Index j) { return this->_values()(i, j); },
        [&](Index i, Index j) { return this->new_dof_B(i, B.permute_ind(j)); });
  }

  /// \brief Return A*config == B*other, store A*config < B*other
  bool operator()(PermuteIterator const &A, PermuteIterator const &B,
                  Eigen::MatrixXd const &other) const {
    _update_A(A, _values());
    _update_B(B, other);
    m_tmp_valid = false;

    return _for_each(
        [&](Index i, Index j) { return this->new_dof_A(i, A.permute_ind(j)); },
        [&](Index i, Index j) { return this->new_dof_B(i, B.permute_ind(j)); });
  }

  /// \brief Returns less than comparison
  ///
  /// - Only valid after call operator returns false
  bool is_less() const { return m_less; }

 private:
  Eigen::MatrixXd const &_values() const { return *m_values_ptr; }

  void _update_A(PermuteIterator const &A,
                 Eigen::MatrixXd const &before) const {
    using clexulator::sublattice_block;
    if (A.factor_group_index() != m_fg_index_A || !m_tmp_valid) {
      for (Index b = 0; b < m_n_sublat; ++b) {
        m_fg_index_A = A.factor_group_index();
        Eigen::MatrixXd const &M = *(A.local_dof_rep(m_key, b).MatrixXd());
        Index dim = M.cols();
        sublattice_block(m_new_dof_A, b, m_n_vol).topRows(dim) =
            M * sublattice_block(before, b, m_n_vol).topRows(dim);
      }
    }
  }

  void _update_B(PermuteIterator const &B,
                 Eigen::MatrixXd const &before) const {
    using clexulator::sublattice_block;
    if (B.factor_group_index() != m_fg_index_B || !m_tmp_valid) {
      for (Index b = 0; b < m_n_sublat; ++b) {
        m_fg_index_B = B.factor_group_index();
        Eigen::MatrixXd const &M = *(B.local_dof_rep(m_key, b).MatrixXd());
        Index dim = M.cols();
        sublattice_block(m_new_dof_B, b, m_n_vol).topRows(dim) =
            M * sublattice_block(before, b, m_n_vol).topRows(dim);
      }
    }
  }

  double new_dof_A(Index i, Index j) const { return m_new_dof_A(i, j); }

  double new_dof_B(Index i, Index j) const { return m_new_dof_B(i, j); }

  template <typename T>
  bool _check(const T &A, const T &B) const {
    // std::cout << "A: " << A << "; B: " << B <<";  ";
    if (A < B - m_tol) {
      m_less = true;
      // std::cout << "Greater!\n";
      return false;
    }
    if (A > B + m_tol) {
      // std::cout << "Less!\n";
      m_less = false;
      return false;
    }
    // std::cout << "Equal!\n";
    return true;
  }

  template <typename F, typename G>
  bool _for_each(F f, G g) const {
    Index i, j;
    for (j = 0; j < _values().cols(); j++) {
      for (i = 0; i < _values().rows(); i++) {
        if (!_check(f(i, j), g(i, j))) {
          return false;
        }
      }
    }
    return true;
  }

  // Points to LocalContinuousConfigDoFValues of ConfigDoF this was constructed
  // with
  Eigen::MatrixXd const *m_values_ptr;

  // DoF type (used to obtain matrix rep)
  DoFKey m_key;

  // Number of sublattices
  Index m_n_sublat;

  // Integer supercell volume
  Index m_n_vol;

  // Tolerance for comparisons
  bool m_tol;

  // Set to false when comparison is made to "other" ConfigDoF, to force update
  // of temporary dof during the next comparison
  mutable bool m_tmp_valid;

  // Store temporary DoFValues under fg operation only:

  mutable Index m_fg_index_A;
  mutable Eigen::MatrixXd m_new_dof_A;

  mutable Index m_fg_index_B;
  mutable Eigen::MatrixXd m_new_dof_B;

  /// Stores (A < B) if A != B
  mutable bool m_less;
};

/// Compare continuous global DoF values
///
/// - Compares global DoF values lexicographically
class Global {
 public:
  Global(Eigen::VectorXd const &_values, DoFKey const &_key, double _tol)
      : m_values_ptr(&_values),
        m_key(_key),
        m_tol(_tol),
        m_tmp_valid(true),
        m_fg_index_A(0),
        m_new_dof_A(*m_values_ptr),
        m_fg_index_B(0),
        m_new_dof_B(*m_values_ptr) {}

  /// \brief Return config == other, store config < other
  bool operator()(Eigen::VectorXd const &other) const {
    return _for_each([&](Index i) { return this->_values(i); },
                     [&](Index i) { return other[i]; });
  }

  /// \brief Return config == B*config, store config < B*config
  bool operator()(PermuteIterator const &B) const {
    _update_B(B, _values());
    m_tmp_valid = true;
    return _for_each([&](Index i) { return this->_values(i); },
                     [&](Index i) { return this->_new_dof_B(i); });
  }

  /// \brief Return A*config == B*config, store A*config < B*config
  bool operator()(PermuteIterator const &A, PermuteIterator const &B) const {
    _update_A(A, _values());
    _update_B(B, _values());
    m_tmp_valid = true;
    return _for_each([&](Index i) { return this->_new_dof_A(i); },
                     [&](Index i) { return this->_new_dof_B(i); });
  }

  /// \brief Return config == B*other, store config < B*other
  bool operator()(PermuteIterator const &B,
                  Eigen::VectorXd const &other) const {
    _update_B(B, other);
    m_tmp_valid = false;
    return _for_each([&](Index i) { return this->_values()[i]; },
                     [&](Index i) { return this->_new_dof_B(i); });
  }

  /// \brief Return A*config == B*other, store A*config < B*other
  bool operator()(PermuteIterator const &A, PermuteIterator const &B,
                  Eigen::VectorXd const &other) const {
    _update_A(A, _values());
    _update_B(B, other);
    m_tmp_valid = false;
    return _for_each([&](Index i) { return this->_new_dof_A(i); },
                     [&](Index i) { return _new_dof_B(i); });
  }

  /// \brief Returns less than comparison
  ///
  /// - Only valid after call operator returns false
  bool is_less() const { return m_less; }

 private:
  void _update_A(PermuteIterator const &A,
                 Eigen::VectorXd const &before) const {
    if (A.factor_group_index() != m_fg_index_A || !m_tmp_valid) {
      m_fg_index_A = A.factor_group_index();
      Eigen::MatrixXd const &M = *(A.global_dof_rep(m_key).MatrixXd());
      m_new_dof_A = M * before;
    }
  }

  void _update_B(PermuteIterator const &B,
                 Eigen::VectorXd const &before) const {
    if (B.factor_group_index() != m_fg_index_B || !m_tmp_valid) {
      m_fg_index_B = B.factor_group_index();
      Eigen::MatrixXd const &M = *(B.global_dof_rep(m_key).MatrixXd());
      m_new_dof_B = M * before;
    }
  }

  Eigen::VectorXd const &_values() const { return *m_values_ptr; }

  double _values(Index i) const { return _values()[i]; }

  double _new_dof_A(Index i) const { return m_new_dof_A[i]; }

  double _new_dof_B(Index i) const { return m_new_dof_B[i]; }

  template <typename T>
  bool _check(const T &A, const T &B) const {
    // std::cout << "A: " << A << "; B: " << B <<";  ";
    if (A < B - m_tol) {
      m_less = true;
      // std::cout << "Greater!\n";
      return false;
    }
    if (A > B + m_tol) {
      // std::cout << "Less!\n";
      m_less = false;
      return false;
    }
    // std::cout << "Equal!\n";
    return true;
  }

  template <typename F, typename G>
  bool _for_each(F f, G g) const {
    Index i;
    for (i = 0; i < _values().size(); i++) {
      if (!_check(f(i), g(i))) {
        return false;
      }
    }
    return true;
  }

  // Points to GlobalContinuousConfigDoFValues of ConfigDoF this was constructed
  // with
  Eigen::VectorXd const *m_values_ptr;

  // DoF type (used to obtain matrix rep)
  DoFKey m_key;

  // Tolerance for comparison
  double m_tol;

  // Set to false when comparison is made to "other" ConfigDoF, to force update
  // of temporary dof during the next comparison
  mutable bool m_tmp_valid;

  // Store temporary DoFValues under fg operation only:

  mutable Index m_fg_index_A;
  mutable Eigen::VectorXd m_new_dof_A;

  mutable Index m_fg_index_B;
  mutable Eigen::VectorXd m_new_dof_B;

  /// Stores (A < B) if A != B
  mutable bool m_less;
};

}  // namespace ConfigDoFIsEquivalent

/** @} */
}  // namespace CASM

#endif
