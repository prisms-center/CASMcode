#include "casm/clex/Clex.hh"

#include "casm/clex/PrimClex.hh"

namespace CASM {


  const PrimClex &Clex::primclex() const {
    return *m_primclex;
  }

  const ClexDescription &Clex::desc() const {
    return m_desc;
  }

  const ClexBasis &Clex::clex_basis() const {

    if(m_clex_basis == nullptr) {
      m_clex_basis = &primclex().clex_basis(m_desc.bset);
    }
    return *m_clex_basis;
  }

  Clexulator &Clex::clexulator(Log &status_log) const {

    if(!m_clexulator.initialized()) {
      m_clexulator = primclex().clexulator(m_desc.bset);
    }
    return m_clexulator;
  }

  const ECIContainer &Clex::eci() const {

    if(m_eci == nullptr) {
      m_eci = &primclex().eci(m_desc);
    }
    return *m_eci;
  }

  /// \brief Compare using descriptions: A.desc() < B.desc()
  bool operator<(const Clex &A, const Clex &B) {
    return A.desc() < B.desc();
  }

}
