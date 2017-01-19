#ifndef CASM_ConfigCompare
#define CASM_ConfigCompare

#include "casm/clex/ConfigIsEquivalent.hh"

namespace CASM {

  /** \ingroup ConfigIsEquivalent
   *
   *  @{
   */

  /// \brief Class for less than comparison of Configurations (with the same Supercell)
  class ConfigCompare {

  public:

    ConfigCompare(const Configuration &_config, double _tol) :
      m_eq(_config, _tol) {}

    /// \brief Check if config < other
    bool operator()(const Configuration &other) const {
      if(m_eq(other)) {
        return false;
      }
      return m_eq.is_less();
    }

    /// \brief Check if config == A*config, store config < A*config
    bool operator()(const PermuteIterator &A) const {
      if(m_eq(A)) {
        return false;
      }
      return m_eq.is_less();
    }

    /// \brief Check if A*config == B*config, store A*config < B*config
    bool operator()(const PermuteIterator &A, const PermuteIterator &B) const {
      if(m_eq(A, B)) {
        return false;
      }
      return m_eq.is_less();
    }

  private:

    ConfigIsEquivalent m_eq;

  };

  /** @} */
}

#endif
