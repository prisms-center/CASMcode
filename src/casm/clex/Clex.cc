#include "casm/clex/Clex.hh"

#include "casm/clex/PrimClex.hh"

namespace CASM {


  const PrimClex &Clex::primclex() const {
    return *m_primclex;
  }

  const ClexDescription &Clex::desc() const {
    return m_desc;
  }

  const SiteOrbitree &Clex::orbitree() const {

    if(m_orbitree == nullptr) {
      m_orbitree = &primclex().orbitree(m_desc);
    }
    return *m_orbitree;
  }

  Clexulator &Clex::clexulator(Log &status_log) const {

    if(!m_clexulator.initialized()) {
      m_clexulator = primclex().clexulator(m_desc);
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
