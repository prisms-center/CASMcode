#ifndef CASM_ConfigIsEquivalent
#define CASM_ConfigIsEquivalent

#include "casm/clex/ConfigDoFIsEquivalent.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"

namespace CASM {

/** \defgroup ConfigIsEquivalent
 *  \ingroup Configuration
 *  \brief Methods for comparing Configurations (with the same Supercell)
 *  @{
 */

/// \brief Class for comparison of Configurations (with the same Supercell)
///
/// - The call operators return the value for equality comparison,
///   and if not equivalent, also store the result for less than comparison
///
class ConfigIsEquivalent {
 public:
  /// Construct with config to be compared against, tolerance for comparison,
  /// and (optional) list of DoFs to compare if _wich_dofs is empty, no dofs
  /// will be compared (default is "all", in which case all DoFs are compared)
  ConfigIsEquivalent(const Configuration &_config, double _tol,
                     std::set<std::string> const &_which_dofs = {"all"})
      : m_config(&_config),
        m_n_sublat(config().configdof().n_sublat()),
        m_all_dofs(_which_dofs.count("all")),
        m_check_occupation(
            (m_all_dofs || _which_dofs.count("occ")) &&
            config().supercell().sym_info().has_occupation_dofs()),
        m_has_aniso_occs(config().supercell().sym_info().has_aniso_occs()),
        m_occupation_ptr(nullptr) {
    ConfigDoF const &configdof = config().configdof();
    clexulator::ConfigDoFValues const &dof_values = configdof.values();

    for (auto const &dof : dof_values.global_dof_values) {
      DoFKey const &key = dof.first;
      Eigen::VectorXd const &values = dof.second;
      if (m_all_dofs || _which_dofs.count(key)) {
        m_global_equivs.emplace(std::piecewise_construct,
                                std::forward_as_tuple(key),
                                std::forward_as_tuple(values, key, _tol));
      }
    }

    if (m_check_occupation) {
      m_occupation_ptr = &dof_values.occupation;
    }

    for (auto const &dof : dof_values.local_dof_values) {
      DoFKey const &key = dof.first;
      Eigen::MatrixXd const &values = dof.second;
      if (m_all_dofs || _which_dofs.count(key)) {
        m_local_equivs.emplace(
            std::piecewise_construct, std::forward_as_tuple(key),
            std::forward_as_tuple(values, key, m_n_sublat, _tol));
      }
    }
  }

  ConfigIsEquivalent(const Configuration &_config,
                     std::set<std::string> const &_which_dofs = {"all"})
      : ConfigIsEquivalent(_config, _config.crystallography_tol(),
                           _which_dofs) {}

  const Configuration &config() const { return *m_config; }

  /// \brief Returns less than comparison
  ///
  /// - Only valid after call operator returns false
  bool is_less() const { return m_less; }

  /// \brief Check if config == other, store config < other
  ///
  /// - Currently assumes that both Configuration have the same Prim, but may
  ///   have different supercells
  bool operator()(const Configuration &other) const {
    if (&config() == &other) {
      return true;
    }

    if (config().shared_prim() != other.shared_prim()) {
      throw std::runtime_error(
          "Error comparing Configuration with ConfigIsEquivalent: "
          "Only Configuration with shared prim may be compared this way.");
    }

    clexulator::ConfigDoFValues const &other_dof_values =
        other.configdof().values();

    if (config().supercell() != other.supercell()) {
      m_less = config().supercell() < other.supercell();
      return false;
    }

    for (auto const &dof_is_equiv_f : m_global_equivs) {
      DoFKey const &key = dof_is_equiv_f.first;
      ConfigDoFIsEquivalent::Global const &f = dof_is_equiv_f.second;
      Eigen::VectorXd const &other_values =
          other_dof_values.global_dof_values.at(key);
      if (!f(other_values)) {
        m_less = f.is_less();
        return false;
      }
    }

    if (!_occupation_is_equivalent(other_dof_values.occupation)) {
      return false;
    }

    for (auto const &dof_is_equiv_f : m_local_equivs) {
      DoFKey const &key = dof_is_equiv_f.first;
      ConfigDoFIsEquivalent::Local const &f = dof_is_equiv_f.second;
      Eigen::MatrixXd const &other_values =
          other_dof_values.local_dof_values.at(key);
      if (!f(other_values)) {
        m_less = f.is_less();
        return false;
      }
    }

    return true;
  }

  /// \brief Check if config == A*config, store config < A*config
  bool operator()(const PermuteIterator &A) const {
    // std::cout << "A.factor_group_index() " << A.factor_group_index() << ";
    // translation_index() " << A.translation_index() << std::endl;

    for (auto const &dof_is_equiv_f : m_global_equivs) {
      ConfigDoFIsEquivalent::Global const &f = dof_is_equiv_f.second;
      if (!f(A)) {
        m_less = f.is_less();
        return false;
      }
    }

    if (!_occupation_is_equivalent(A)) {
      return false;
    }

    for (auto const &dof_is_equiv_f : m_local_equivs) {
      ConfigDoFIsEquivalent::Local const &f = dof_is_equiv_f.second;
      if (!f(A)) {
        m_less = f.is_less();
        return false;
      }
    }

    return true;
  }

  /// \brief Check if A*config == B*config, store A*config < B*config
  bool operator()(const PermuteIterator &A, const PermuteIterator &B) const {
    if (A.factor_group_index() != B.factor_group_index()) {
      for (auto const &dof_is_equiv_f : m_global_equivs) {
        ConfigDoFIsEquivalent::Global const &f = dof_is_equiv_f.second;
        if (!f(A, B)) {
          m_less = f.is_less();
          return false;
        }
      }
    }

    if (!_occupation_is_equivalent(A, B)) {
      return false;
    }

    for (auto const &dof_is_equiv_f : m_local_equivs) {
      ConfigDoFIsEquivalent::Local const &f = dof_is_equiv_f.second;
      if (!f(A, B)) {
        m_less = f.is_less();
        return false;
      }
    }

    return true;
  }

  /// \brief Check if config == A*other, store config < A*other
  bool operator()(const PermuteIterator &A, const Configuration &other) const {
    clexulator::ConfigDoFValues const &other_dof_values =
        other.configdof().values();

    for (auto const &dof_is_equiv_f : m_global_equivs) {
      DoFKey const &key = dof_is_equiv_f.first;
      ConfigDoFIsEquivalent::Global const &f = dof_is_equiv_f.second;
      Eigen::VectorXd const &other_values =
          other_dof_values.global_dof_values.at(key);
      if (!f(A, other_values)) {
        m_less = f.is_less();
        return false;
      }
    }

    if (!_occupation_is_equivalent(A, other_dof_values.occupation)) {
      return false;
    }

    for (auto const &dof_is_equiv_f : m_local_equivs) {
      DoFKey const &key = dof_is_equiv_f.first;
      ConfigDoFIsEquivalent::Local const &f = dof_is_equiv_f.second;
      Eigen::MatrixXd const &other_values =
          other_dof_values.local_dof_values.at(key);
      if (!f(A, other_values)) {
        m_less = f.is_less();
        return false;
      }
    }

    return true;
  }

  /// \brief Check if A*config == B*other, store A*config < B*other
  bool operator()(const PermuteIterator &A, const PermuteIterator &B,
                  const Configuration &other) const {
    clexulator::ConfigDoFValues const &other_dof_values =
        other.configdof().values();

    for (auto const &dof_is_equiv_f : m_global_equivs) {
      DoFKey const &key = dof_is_equiv_f.first;
      ConfigDoFIsEquivalent::Global const &f = dof_is_equiv_f.second;
      Eigen::VectorXd const &other_values =
          other_dof_values.global_dof_values.at(key);
      if (!f(A, B, other_values)) {
        m_less = f.is_less();
        return false;
      }
    }

    if (!_occupation_is_equivalent(A, B, other_dof_values.occupation)) {
      return false;
    }

    for (auto const &dof_is_equiv_f : m_local_equivs) {
      DoFKey const &key = dof_is_equiv_f.first;
      ConfigDoFIsEquivalent::Local const &f = dof_is_equiv_f.second;
      Eigen::MatrixXd const &other_values =
          other_dof_values.local_dof_values.at(key);
      if (!f(A, B, other_values)) {
        m_less = f.is_less();
        return false;
      }
    }

    return true;
  }

 private:
  template <typename... Args>
  bool _occupation_is_equivalent(Args &&...args) const;

  const Configuration *m_config;
  Index m_n_sublat;
  bool m_all_dofs;
  bool m_check_occupation;
  bool m_has_aniso_occs;
  Eigen::VectorXi const *m_occupation_ptr;
  std::map<DoFKey, ConfigDoFIsEquivalent::Global> m_global_equivs;
  std::map<DoFKey, ConfigDoFIsEquivalent::Local> m_local_equivs;
  mutable bool m_less;
};

template <typename... Args>
bool ConfigIsEquivalent::_occupation_is_equivalent(Args &&...args) const {
  if (m_check_occupation) {
    if (m_has_aniso_occs) {
      ConfigDoFIsEquivalent::AnisoOccupation f(*m_occupation_ptr, m_n_sublat);
      if (!f(std::forward<Args>(args)...)) {
        m_less = f.is_less();
        return false;
      }
    } else {
      ConfigDoFIsEquivalent::Occupation f(*m_occupation_ptr);
      if (!f(std::forward<Args>(args)...)) {
        m_less = f.is_less();
        return false;
      }
    }
  }
  return true;
}

/** @} */
}  // namespace CASM

#endif
