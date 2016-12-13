#ifndef CASM_ConfigDoFCompare
#define CASM_ConfigDoFCompare

#include "casm/clex/ConfigDoFIsEquivalent.hh"

namespace CASM {

  /** \ingroup ConfigIsEquivalent
   *
   *  @{
   */

  /// \brief Wrapper class for generic less than comparison of ConfigDoF
  ///
  /// - Wraps a functor derived from ConfigDoFIsEquivalentBase that is specialized for
  ///   comparison of a particular type of DoF
  class ConfigDoFCompare {

  public:

    /// \brief Construct a ConfigDoFCompare object for a particular DoF type
    ///
    /// Easiest construction is probably using 'make_dof_compare'.
    ///
    /// Example:
    /// \code
    /// ConfigDoFCompare strain_compare = make_dof_compare<DoFIsEquivalent::Strain>(my_configdof or my_config);
    /// ConfigDoFCompare occ_compare = make_dof_compare<DoFIsEquivalent::Occupation>(my_configdof or my_config);
    /// ConfigDoFCompare disp_compare = make_dof_compare<DoFIsEquivalent::Displacement>(my_configdof or my_config);
    /// \endcode
    template<typename ConfigDoFIsEquivalentType>
    ConfigDoFCompare(std::unique_ptr<ConfigDoFIsEquivalentType> f) :
      m_f(f) {}


    /// \brief Return config < other
    bool operator()(const Configuration &other) const {
      return (*this)(other.configdof());
    }

    /// \brief Return config < other
    bool operator()(const ConfigDoF &other) const {
      if((*m_f)(other)) {
        return false;
      }
      return m_f->is_less();
    }

    /// \brief Return config < A*config
    bool operator()(const PermuteIterator &A) const {
      if((*m_f)(A)) {
        return false;
      }
      return m_f->is_less();
    }

    /// \brief Return A*config < B*config
    bool operator()(const PermuteIterator &A, const PermuteIterator &B) const {
      if((*m_f)(A, B)) {
        return false;
      }
      return m_f->is_less();
    }

  private:
    notstd::cloneable_ptr<DoFIsEquivalent::ConfigDoFIsEquivalentBase> m_f;

  };

  /// Factory function to make ConfigDoFCompare
  ///
  /// \relates ConfigDoFCompare
  template<typename ConfigDoFIsEquivalentType, typename ...Args>
  ConfigDoFCompare make_dof_compare(Args &&...args) {
    return ConfigDoFCompare(notstd::make_unique<ConfigDoFIsEquivalentType>(std::forward<Args>(args)...));
  }

  /** @} */
}

#endif
