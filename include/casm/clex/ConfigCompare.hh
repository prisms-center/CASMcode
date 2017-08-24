#ifndef CASM_ConfigCompare
#define CASM_ConfigCompare

#include <utility>

namespace CASM {

  class PermuteIterator;

  /** \ingroup ConfigIsEquivalent
   *
   *  @{
   */

  /// \brief Class for less than comparison of Configurations implemented via
  ///        a ConfigTypeIsEqual class that also stores the less than result
  template<typename ConfigType, typename IsEqualImpl>
  class GenericConfigCompare {

  public:

    explicit GenericConfigCompare(const IsEqualImpl &_eq) :
      m_eq(_eq) {}

    template<typename... Args>
    bool operator()(Args &&... args) const {
      if(m_eq(std::forward<Args>(args)...)) {
        return false;
      }
      return m_eq.is_less();
    }

    /*
    /// \brief Return config < other (may have different Supercell)
    bool operator()(const ConfigType &other) const {
      if(&m_eq.config().supercell() != &other.supercell()) {
        if(m_eq.config().supercell() != other.supercell()) {
          return m_eq.config().supercell() < other.supercell();
        }
      }
      if(m_eq(other)) {
        return false;
      }
      return m_eq.is_less();
    }

    /// \brief Return config < A*config
    bool operator()(const PermuteIterator &A) const {
      if(m_eq(A)) {
        return false;
      }
      return m_eq.is_less();
    }

    /// \brief Return A*config < B*config
    bool operator()(const PermuteIterator &A, const PermuteIterator &B) const {
      if(m_eq(A, B)) {
        return false;
      }
      return m_eq.is_less();
    }

    /// \brief Return config < A*other
    bool operator()(const PermuteIterator &A, const ConfigType& other) const {
      if(m_eq(A, other)) {
        return false;
      }
      return m_eq.is_less();
    }

    /// \brief Return A*config < B*other
    bool operator()(const PermuteIterator &A, const PermuteIterator &B, const ConfigType& other) const {
      if(m_eq(A, B, other)) {
        return false;
      }
      return m_eq.is_less();
    }
    */

    const IsEqualImpl &base() const {
      return m_eq;
    }

  private:

    IsEqualImpl m_eq;

  };

  /** @} */
}

#endif
