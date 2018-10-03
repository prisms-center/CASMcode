#ifndef CASM_ConfigIsEquivalent
#define CASM_ConfigIsEquivalent

#include "casm/clex/ConfigDoFCompare.hh"
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
    using EquivPtr = notstd::cloneable_ptr<ConfigDoFIsEquivalent::Base>;

    ConfigIsEquivalent(const Configuration &_config, double _tol) :
      m_config(&_config) {

      for(auto const &dof : config().configdof().global_dofs())
        m_global_equivs.push_back(notstd::make_cloneable<ConfigDoFIsEquivalent::Global>(_config, dof.first, _tol));

      if(config().configdof().has_occupation()) {
        m_local_equivs.push_back(notstd::make_cloneable<ConfigDoFIsEquivalent::Occupation>(_config.configdof()));
      }

      for(auto const &dof : config().configdof().local_dofs())
        m_local_equivs.push_back(notstd::make_cloneable<ConfigDoFIsEquivalent::Local>(_config, dof.first, _tol));

    }

    ConfigIsEquivalent(const Configuration &_config) :
      ConfigIsEquivalent(_config, _config.crystallography_tol()) {}

    const Configuration &config() const {
      return *m_config;
    }

    std::vector<EquivPtr> const &local_equivs() const {
      return m_local_equivs;
    }

    std::vector<EquivPtr> const &global_equivs() const {
      return m_global_equivs;
    }


    /// \brief Returns less than comparison
    ///
    /// - Only valid after call operator returns false
    bool is_less() const {
      return m_less;
    }

    /// \brief Check if config == other, store config < other
    ///
    /// - Currently assumes that both Configuration have the same DoF types
    bool operator()(const Configuration &other) const {
      if(&config() == &other) {
        return true;
      }

      if(config().supercell() != other.supercell()) {
        m_less = config().supercell() < other.supercell();
        return false;
      }

      for(const auto &eq : global_equivs()) {
        // check if config == other, for this global DoF type
        if(!(*eq)(other.configdof())) {
          m_less = eq->is_less();
          return false;
        }
      }

      for(const auto &eq : local_equivs()) {
        // check if config == other, for this local DoF type
        if(!(*eq)(other.configdof())) {
          m_less = eq->is_less();
          return false;
        }
      }

      return true;
    }

    /// \brief Check if config == A*config, store config < A*config
    bool operator()(const PermuteIterator &A) const {
      for(const auto &eq : global_equivs()) {
        // check if config == A*config, for this global DoF type
        if(!(*eq)(A)) {
          m_less = eq->is_less();
          return false;
        }
      }

      for(const auto &eq : local_equivs()) {
        // check if config == A*config, for this local DoF type
        if(!(*eq)(A)) {
          m_less = eq->is_less();
          return false;
        }
      }

      return true;
    }

    /// \brief Check if A*config == B*config, store A*config < B*config
    bool operator()(const PermuteIterator &A, const PermuteIterator &B) const {

      if(A.factor_group_index() != B.factor_group_index()) {
        for(const auto &eq : global_equivs()) {
          // check if config == other, for this global DoF type
          if(!(*eq)(A, B)) {
            m_less = eq->is_less();
            return false;
          }
        }
      }

      for(const auto &eq : local_equivs()) {
        // check if config == other, for this local DoF type
        if(!(*eq)(A, B)) {
          m_less = eq->is_less();
          return false;
        }
      }

      return true;
    }

    /// \brief Check if config == A*other, store config < A*other
    bool operator()(const PermuteIterator &A, const Configuration &other) const {

      for(const auto &eq : global_equivs()) {
        // check if config == other, for this global DoF type
        if(!(*eq)(A, other.configdof())) {
          m_less = eq->is_less();
          return false;
        }
      }

      for(const auto &eq : local_equivs()) {
        // check if config == other, for this local DoF type
        if(!(*eq)(A, other.configdof())) {
          m_less = eq->is_less();
          return false;
        }
      }

      return true;
    }

    /// \brief Check if A*config == B*other, store A*config < B*other
    bool operator()(const PermuteIterator &A, const PermuteIterator &B, const Configuration &other) const {

      for(const auto &eq : global_equivs()) {
        // check if config == other, for this global DoF type
        if(!(*eq)(A, B, other.configdof())) {
          m_less = eq->is_less();
          return false;
        }
      }

      for(const auto &eq : local_equivs()) {
        // check if config == other, for this local DoF type
        if(!(*eq)(A, B, other.configdof())) {
          m_less = eq->is_less();
          return false;
        }
      }

      return true;
    }

  private:

    const Configuration *m_config;
    std::vector<EquivPtr> m_global_equivs;
    std::vector<EquivPtr> m_local_equivs;
    mutable bool m_less;
  };

  /** @} */
}

#endif
