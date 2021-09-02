#ifndef CASM_ConfigDoFIsEquivalent
#define CASM_ConfigDoFIsEquivalent

#include "casm/casm_io/Log.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clexulator/ConfigDoFValuesTools.hh"

namespace CASM {

/** \ingroup ConfigIsEquivalent
 *  @{
 */

/// Namespace containing DoF comparison functors
namespace ConfigDoFIsEquivalent {

/// \brief Base class for functors that compare ConfigDoF
///
/// - Specialized for Occupation or Float DoF, and then for type
///   (strain, occupation, displacement etc.)
/// - The derived call operators return the value for equality comparison,
///   and if not equivalent, also store the result for less than comparison
class Base {
 public:
  virtual ~Base() {}

  /// \brief Returns less than comparison
  ///
  /// - Only valid after call operator returns false
  bool is_less() const { return m_less; }

  /// \brief Return config == other
  bool operator()(Configuration const &other) const {
    return (*this)(other.configdof());
  }

  /// \brief Return config == other
  virtual bool operator()(ConfigDoF const &other) const = 0;

  /// \brief Return config == A*config
  virtual bool operator()(PermuteIterator const &A) const = 0;

  /// \brief Return A*config == B*config
  virtual bool operator()(PermuteIterator const &A,
                          PermuteIterator const &B) const = 0;

  /// \brief Return config == A*other
  virtual bool operator()(PermuteIterator const &A,
                          ConfigDoF const &other) const = 0;

  /// \brief Return A*config == B*other
  virtual bool operator()(PermuteIterator const &A, PermuteIterator const &B,
                          ConfigDoF const &other) const = 0;

  std::unique_ptr<Base> clone() const {
    return std::unique_ptr<Base>(this->_clone());
  }

 protected:
  mutable bool m_less;

 private:
  virtual Base *_clone() const = 0;
};

/// Compare isotropic occupation values
///
/// - The protected '_check' method provides for both checking equality and if
///   not equivalent, storing the 'less than' result
class Occupation : public Base {
 public:
  Occupation(ConfigDoF const &_configdof) : m_configdof_ptr(&_configdof) {}

  /// \brief Return config == other, store config < other
  bool operator()(ConfigDoF const &other) const override {
    return _for_each([&](Index i) { return this->configdof().occ(i); },
                     [&](Index i) { return other.occ(i); });
  }

  /// \brief Return config == A*config, store config < A*config
  bool operator()(PermuteIterator const &A) const override {
    return _for_each(
        [&](Index i) { return this->configdof().occ(i); },
        [&](Index i) { return this->configdof().occ(A.permute_ind(i)); });
  }

  /// \brief Return A*config == B*config, store A*config < B*config
  bool operator()(PermuteIterator const &A,
                  PermuteIterator const &B) const override {
    return _for_each(
        [&](Index i) { return this->configdof().occ(A.permute_ind(i)); },
        [&](Index i) { return this->configdof().occ(B.permute_ind(i)); });
  }

  /// \brief Return config == A*other, store config < A*other
  bool operator()(PermuteIterator const &A,
                  ConfigDoF const &other) const override {
    return _for_each([&](Index i) { return this->configdof().occ(i); },
                     [&](Index i) { return other.occ(A.permute_ind(i)); });
  }

  /// \brief Return A*config == B*other, store A*config < B*other
  bool operator()(PermuteIterator const &A, PermuteIterator const &B,
                  ConfigDoF const &other) const override {
    return _for_each(
        [&](Index i) { return this->configdof().occ(A.permute_ind(i)); },
        [&](Index i) { return other.occ(B.permute_ind(i)); });
  }

 protected:
  ConfigDoF const &configdof() const { return *m_configdof_ptr; }

  template <typename F, typename G>
  bool _for_each(F f, G g) const {
    Index i;
    for (i = 0; i < configdof().size(); i++) {
      if (!_check(f(i), g(i))) {
        return false;
      }
    }
    // std::cout << "Equal!\n";
    return true;
  }

 private:
  Occupation *_clone() const override { return new Occupation(*this); }

  template <typename T>
  bool _check(const T &A, const T &B) const {
    if (A == B) {
      return true;
    }
    m_less = (A < B);
    return false;
  }

  ConfigDoF const *m_configdof_ptr;
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
class AnisoOccupation : public Base {
 public:
  AnisoOccupation(ConfigDoF const &_configdof)
      : m_configdof_ptr(&_configdof),
        m_tmp_valid(true),
        m_fg_index_A(0),
        m_new_occ_A(m_configdof_ptr->occupation()),
        m_fg_index_B(0),
        m_new_occ_B(m_configdof_ptr->occupation()) {}

  /// \brief Return config == other, store config < other
  bool operator()(ConfigDoF const &other) const override {
    return _for_each([&](Index i) { return this->configdof().occ(i); },
                     [&](Index i) { return other.occ(i); });
  }

  /// \brief Return config == B*config, store config < B*config
  bool operator()(PermuteIterator const &B) const override {
    _update_B(B, configdof());
    m_tmp_valid = true;

    return _for_each(
        [&](Index i) { return this->configdof().occ(i); },
        [&](Index i) { return this->m_new_occ_B[B.permute_ind(i)]; });
  }

  /// \brief Return A*config == B*config, store A*config < B*config
  bool operator()(PermuteIterator const &A,
                  PermuteIterator const &B) const override {
    _update_A(A, configdof());
    _update_B(B, configdof());
    m_tmp_valid = true;

    return _for_each(
        [&](Index i) { return this->m_new_occ_A[A.permute_ind(i)]; },
        [&](Index i) { return this->m_new_occ_B[B.permute_ind(i)]; });
  }

  /// \brief Return config == B*other, store config < B*other
  bool operator()(PermuteIterator const &B,
                  ConfigDoF const &other) const override {
    _update_B(B, other);
    m_tmp_valid = false;

    return _for_each(
        [&](Index i) { return this->m_configdof_ptr->occ(i); },
        [&](Index i) { return this->m_new_occ_B[B.permute_ind(i)]; });
  }

  /// \brief Return A*config == B*other, store A*config < B*other
  bool operator()(PermuteIterator const &A, PermuteIterator const &B,
                  ConfigDoF const &other) const override {
    _update_A(A, configdof());
    _update_B(B, other);
    m_tmp_valid = false;

    return _for_each(
        [&](Index i) { return this->m_new_occ_A[A.permute_ind(i)]; },
        [&](Index i) { return this->m_new_occ_B[B.permute_ind(i)]; });
  }

 protected:
  ConfigDoF const &configdof() const { return *m_configdof_ptr; }

  template <typename F, typename G>
  bool _for_each(F f, G g) const {
    Index i;
    for (i = 0; i < configdof().size(); i++) {
      if (!_check(f(i), g(i))) {
        return false;
      }
    }
    // std::cout << "Equal!\n";
    return true;
  }

  void _update_A(PermuteIterator const &A, ConfigDoF const &before) const {
    if (A.factor_group_index() != m_fg_index_A || !m_tmp_valid) {
      m_fg_index_A = A.factor_group_index();
      Index l = 0;
      for (Index b = 0; b < m_configdof_ptr->n_sublat(); ++b) {
        for (Index n = 0; n < m_configdof_ptr->n_vol(); ++n, ++l) {
          m_new_occ_A[l] = (*(A.occ_rep(b).permutation()))[before.occ(l)];
        }
      }
    }
  }

  void _update_B(PermuteIterator const &B, ConfigDoF const &before) const {
    if (B.factor_group_index() != m_fg_index_B || !m_tmp_valid) {
      m_fg_index_B = B.factor_group_index();
      Index l = 0;
      for (Index b = 0; b < m_configdof_ptr->n_sublat(); ++b) {
        for (Index n = 0; n < m_configdof_ptr->n_vol(); ++n, ++l) {
          m_new_occ_B[l] = (*(B.occ_rep(b).permutation()))[before.occ(l)];
        }
      }
    }
  }

 private:
  AnisoOccupation *_clone() const override {
    return new AnisoOccupation(*this);
  }

  template <typename T>
  bool _check(const T &A, const T &B) const {
    if (A == B) {
      return true;
    }
    m_less = (A < B);
    return false;
  }

  // Points to ConfigDoF this was constructed with
  ConfigDoF const *m_configdof_ptr;

  // Set to false when comparison is made to "other" ConfigDoF, to force update
  // of temporary dof during the next comparison
  mutable bool m_tmp_valid;

  // Store temporary DoFValues under fg operation only:

  mutable Index m_fg_index_A;
  mutable Eigen::VectorXi m_new_occ_A;

  mutable Index m_fg_index_B;
  mutable Eigen::VectorXi m_new_occ_B;
};

/// \brief Abstract base class specialization of Base for
/// floating point DoF types
///
/// - The protected '_check' method provides for both checking equality and if
///   not equivalent, storing the 'less than' result
class Float : public Base {
 public:
  Float(double _tol, DoFKey const &_key) : m_tol(_tol), m_key(_key) {}

  DoFKey const &key() const { return m_key; }

 protected:
  template <typename T>
  bool _check(const T &A, const T &B) const {
    // std::cout << "A: " << A << "; B: " << B <<";  ";
    if (A < B - tol()) {
      m_less = true;
      // std::cout << "Greater!\n";
      return false;
    }
    if (A > B + tol()) {
      // std::cout << "Less!\n";
      m_less = false;
      return false;
    }
    // std::cout << "Equal!\n";
    return true;
  }

 private:
  double tol() const { return m_tol; }

  const double m_tol;

  const DoFKey m_key;
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
class Local : public Float {
 public:
  Local(ConfigDoF const &_configdof, DoFKey const &_key, double _tol)
      : Float(_tol, _key),
        m_n_sublat(_configdof.n_sublat()),
        m_n_vol(_configdof.n_vol()),
        m_values_ptr(&(_configdof.local_dof(_key).values())),
        m_tmp_valid(true),
        m_zeros(*m_values_ptr),
        m_fg_index_A(0),
        m_new_dof_A(*m_values_ptr),
        m_fg_index_B(0),
        m_new_dof_B(*m_values_ptr) {
    m_zeros.setZero();
  }

  Local(Configuration const &_config, DoFKey const &_key, double _tol)
      : Local(_config.configdof(), _key, _tol) {}

  /// \brief Return config == other, store config < other
  bool operator()(ConfigDoF const &other) const override {
    Eigen::MatrixXd const *other_ptr;
    if (other.has_local_dof(key())) {
      other_ptr = &(other.local_dof(key()).values());
    } else {
      other_ptr = &m_zeros;
    }

    return _for_each(
        [&](Index i, Index j) { return this->_values()(i, j); },
        [&](Index i, Index j) { return (*other_ptr)(i, j); });
  }

  /// \brief Return config == B*config, store config < B*config
  bool operator()(PermuteIterator const &B) const override {
    _update_B(B, _values());
    m_tmp_valid = true;

    return _for_each(
        [&](Index i, Index j) { return this->_values()(i, j); },
        [&](Index i, Index j) { return this->new_dof_B(i, B.permute_ind(j)); });
  }

  /// \brief Return A*config == B*config, store A*config < B*config
  bool operator()(PermuteIterator const &A,
                  PermuteIterator const &B) const override {
    _update_A(A, _values());
    _update_B(B, _values());
    m_tmp_valid = true;
    return _for_each(
        [&](Index i, Index j) { return this->new_dof_A(i, A.permute_ind(j)); },
        [&](Index i, Index j) { return this->new_dof_B(i, B.permute_ind(j)); });
  }

  /// \brief Return config == B*other, store config < B*other
  bool operator()(PermuteIterator const &B,
                  ConfigDoF const &other) const override {
    Eigen::MatrixXd const *other_ptr;
    if (other.has_local_dof(key())) {
      other_ptr = &(other.local_dof(key()).values());
    } else {
      other_ptr = &m_zeros;
    }

    _update_B(B, *other_ptr);
    m_tmp_valid = false;

    return _for_each(
        [&](Index i, Index j) { return this->_values()(i, j); },
        [&](Index i, Index j) { return this->new_dof_B(i, B.permute_ind(j)); });
  }

  /// \brief Return A*config == B*other, store A*config < B*other
  bool operator()(PermuteIterator const &A, PermuteIterator const &B,
                  ConfigDoF const &other) const override {
    Eigen::MatrixXd const *other_ptr;
    if (other.has_local_dof(key())) {
      other_ptr = &(other.local_dof(key()).values());
    } else {
      other_ptr = &m_zeros;
    }
    _update_A(A, _values());
    _update_B(B, *other_ptr);
    m_tmp_valid = false;

    return _for_each(
        [&](Index i, Index j) { return this->new_dof_A(i, A.permute_ind(j)); },
        [&](Index i, Index j) { return this->new_dof_B(i, B.permute_ind(j)); });
  }

 private:
  Eigen::MatrixXd const &_values() const {
    return *m_values_ptr;
  }

  void _update_A(PermuteIterator const &A,
                 Eigen::MatrixXd const &before) const {
    using clexulator::sublattice_block;
    if (A.factor_group_index() != m_fg_index_A || !m_tmp_valid) {
      for (Index b = 0; b < m_n_sublat; ++b) {
        m_fg_index_A = A.factor_group_index();
        Eigen::MatrixXd const &M = *(A.local_dof_rep(key(), b).MatrixXd());
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
        Eigen::MatrixXd const &M = *(B.local_dof_rep(key(), b).MatrixXd());
        Index dim = M.cols();
        sublattice_block(m_new_dof_B, b, m_n_vol).topRows(dim) =
            M * sublattice_block(before, b, m_n_vol).topRows(dim);
      }
    }
  }

  double new_dof_A(Index i, Index j) const {
    return m_new_dof_A(i, j);
  }

  double new_dof_B(Index i, Index j) const {
    return m_new_dof_B(i, j);
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

  Base *_clone() const override { return new Local(*this); }

  // Number of sublattices
  Index m_n_sublat;

  // Integer supercell volume
  Index m_n_vol;

  // Points to LocalContinuousConfigDoFValues of ConfigDoF this was constructed
  // with
  Eigen::MatrixXd const *m_values_ptr;

  // Set to false when comparison is made to "other" ConfigDoF, to force update
  // of temporary dof during the next comparison
  mutable bool m_tmp_valid;

  // Used when "other" ConfigDoF does not have local DoF==this->key()
  mutable Eigen::MatrixXd m_zeros;

  // Store temporary DoFValues under fg operation only:

  mutable Index m_fg_index_A;
  mutable Eigen::MatrixXd m_new_dof_A;

  mutable Index m_fg_index_B;
  mutable Eigen::MatrixXd m_new_dof_B;
};

/// Compare continuous global DoF values
///
/// - Compares global DoF values lexicographically
class Global : public Float {
 public:
  Global(ConfigDoF const &_configdof, DoFKey const &_key, double _tol)
      : Float(_tol, _key),
        m_values_ptr(&_configdof.global_dof(_key).values()),
        m_tmp_valid(true),
        m_zeros(*m_values_ptr),
        m_fg_index_A(0),
        m_new_dof_A(*m_values_ptr),
        m_fg_index_B(true),
        m_new_dof_B(*m_values_ptr) {
    m_zeros.setZero();
  }

  Global(Configuration const &_config, DoFKey const &_key, double _tol)
      : Global(_config.configdof(), _key, _tol) {}

  /// \brief Return config == other, store config < other
  bool operator()(ConfigDoF const &other) const override {
    Eigen::VectorXd const *other_ptr;
    if (other.has_global_dof(key())) {
      other_ptr = &(other.global_dof(key()).values());
    } else {
      other_ptr = &m_zeros;
    }

    return _for_each([&](Index i) { return this->_values(i); },
                     [&](Index i) { return (*other_ptr)[i]; });
  }

  /// \brief Return config == B*config, store config < B*config
  bool operator()(PermuteIterator const &B) const override {
    _update_B(B, _values());
    m_tmp_valid = true;
    return _for_each([&](Index i) { return this->_values(i); },
                     [&](Index i) { return this->_new_dof_B(i); });
  }

  /// \brief Return A*config == B*config, store A*config < B*config
  bool operator()(PermuteIterator const &A,
                  PermuteIterator const &B) const override {
    _update_A(A, _values());
    _update_B(B, _values());
    m_tmp_valid = true;
    return _for_each([&](Index i) { return this->_new_dof_A(i); },
                     [&](Index i) { return this->_new_dof_B(i); });
  }

  /// \brief Return config == B*other, store config < B*other
  bool operator()(PermuteIterator const &B,
                  ConfigDoF const &other) const override {
    Eigen::VectorXd const *other_ptr;
    if (other.has_global_dof(key())) {
      other_ptr = &(other.global_dof(key()).values());
    } else {
      other_ptr = &m_zeros;
    }
    _update_B(B, *other_ptr);
    m_tmp_valid = false;
    return _for_each([&](Index i) { return this->_values()[i]; },
                     [&](Index i) { return this->_new_dof_B(i); });
  }

  /// \brief Return A*config == B*other, store A*config < B*other
  bool operator()(PermuteIterator const &A, PermuteIterator const &B,
                  ConfigDoF const &other) const override {
    Eigen::VectorXd const *other_ptr;
    if (other.has_global_dof(key())) {
      other_ptr = &(other.global_dof(key()).values());
    } else {
      other_ptr = &m_zeros;
    }
    _update_A(A, _values());
    _update_B(B, *other_ptr);
    m_tmp_valid = false;
    return _for_each([&](Index i) { return this->_new_dof_A(i); },
                     [&](Index i) { return _new_dof_B(i); });
  }

 private:
  Base *_clone() const override { return new Global(*this); }

  void _update_A(PermuteIterator const &A,
                 Eigen::VectorXd const &before) const {
    if (A.factor_group_index() != m_fg_index_A || !m_tmp_valid) {
      m_fg_index_A = A.factor_group_index();
      Eigen::MatrixXd const &M = *(A.global_dof_rep(key()).MatrixXd());
      m_new_dof_A = M * before;
    }
  }

  void _update_B(PermuteIterator const &B,
                 Eigen::VectorXd const &before) const {
    if (B.factor_group_index() != m_fg_index_B || !m_tmp_valid) {
      m_fg_index_B = B.factor_group_index();
      Eigen::MatrixXd const &M = *(B.global_dof_rep(key()).MatrixXd());
      m_new_dof_B = M * before;
    }
  }

  Eigen::VectorXd const &_values() const {
    return *m_values_ptr;
  }

  double _values(Index i) const { return _values()[i]; }

  double _new_dof_A(Index i) const { return m_new_dof_A[i]; }

  double _new_dof_B(Index i) const { return m_new_dof_B[i]; }

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

  // Set to false when comparison is made to "other" ConfigDoF, to force update
  // of temporary dof during the next comparison
  mutable bool m_tmp_valid;

  // Used when "other" ConfigDoF does not have local DoF==this->key()
  mutable Eigen::VectorXd m_zeros;

  // Store temporary DoFValues under fg operation only:

  mutable Index m_fg_index_A;
  mutable Eigen::VectorXd m_new_dof_A;

  mutable Index m_fg_index_B;
  mutable Eigen::VectorXd m_new_dof_B;
};

}  // namespace ConfigDoFIsEquivalent

/// Factory function to make ConfigDoFIsEquivalent
///
/// \relates ConfigDoFIsEquivalent
// template<typename ConfigDoFIsEquivalentType, typename ...Args>
// ConfigDoFIsEquivalent make_dof_is_equivalent(Args &&...args) {
// return
// ConfigDoFIsEquivalent(notstd::make_unique<ConfigDoFIsEquivalentType>(std::forward<Args>(args)...));
//}

/** @} */
}  // namespace CASM

#endif
