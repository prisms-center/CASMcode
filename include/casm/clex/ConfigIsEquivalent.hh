#ifndef CASM_ConfigIsEquivalent
#define CASM_ConfigIsEquivalent

#include "casm/clex/ConfigDoFCompare.hh"


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

    typedef std::vector<ConfigDoFIsEquivalent> eq_container;

    ConfigIsEquivalent(const Configuration &_config, double _tol) :
      m_config(&_config) {

      if(config().has_deformation()) {
        m_global_eq.push_back(make_dof_is_equivalent<DoFIsEquivalent::Strain>(_config, _tol));
      }

      if(config().has_occupation()) {
        m_site_eq.push_back(make_dof_is_equivalent<DoFIsEquivalent::Occupation>(_config));
      }

      if(config().has_displacement()) {
        m_site_eq.push_back(make_dof_is_equivalent<DoFIsEquivalent::Displacement>(_config, _tol));
      }
    }

    const Configuration &config() const {
      return *m_config;
    }

    eq_container &global_eq() {
      return m_global_eq;
    }

    const eq_container &global_eq() const {
      return m_global_eq;
    }

    eq_container &site_eq() {
      return m_site_eq;
    }

    const eq_container &site_eq() const {
      return m_site_eq;
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

      for(const auto &g : global_eq()) {
        // check if config == other, for this global DoF type
        if(!g(other.configdof())) {
          m_less = g.is_less();
          return false;
        }
      }

      for(const auto &s : site_eq()) {
        // check if config == other, for this site DoF type
        if(!s(other.configdof())) {
          m_less = s.is_less();
          return false;
        }
      }

      return true;
    }

    /// \brief Check if config == A*config, store config < A*config
    bool operator()(const PermuteIterator &A) const {
      for(const auto &g : global_eq()) {
        // check if config == other, for this global DoF type
        if(!g(A)) {
          m_less = g.is_less();
          return false;
        }
      }

      for(const auto &s : site_eq()) {
        // check if config == other, for this site DoF type
        if(!s(A)) {
          m_less = s.is_less();
          return false;
        }
      }

      return true;
    }

    /// \brief Check if A*config == B*config, store A*config < B*config
    bool operator()(const PermuteIterator &A, const PermuteIterator &B) const {

      if(A.factor_group_index() != B.factor_group_index()) {
        for(const auto &g : global_eq()) {
          // check if config == other, for this global DoF type
          if(!g(A, B)) {
            m_less = g.is_less();
            return false;
          }
        }
      }

      for(const auto &s : site_eq()) {
        // check if config == other, for this site DoF type
        if(!s(A, B)) {
          m_less = s.is_less();
          return false;
        }
      }

      return true;
    }

  private:

    const Configuration *m_config;
    eq_container m_global_eq;
    eq_container m_site_eq;
    mutable bool m_less;
  };

  /** @} */
}

#endif
