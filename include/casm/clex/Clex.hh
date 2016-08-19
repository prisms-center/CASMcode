#ifndef CASM_Clex
#define CASM_Clex
#include <cstddef>

#include "casm/CASM_global_definitions.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/app/ProjectSettings.hh"

namespace CASM {

  class PrimClex;
  class ECIContainer;
  class ClexBasis;

  /// \brief Data structure used for cluster expansions
  class Clex {

  public:

    Clex() {}

    Clex(const PrimClex &_primclex) :
      m_primclex(&_primclex), m_clex_basis(nullptr), m_eci(nullptr) {}

    Clex(const PrimClex &_primclex, const ClexDescription &_desc) :
      m_primclex(&_primclex), m_desc(_desc), m_clex_basis(nullptr), m_eci(nullptr) {}


    const PrimClex &primclex() const;

    const ClexDescription &desc() const;

    const ClexBasis &clex_basis() const;

    Clexulator &clexulator(Log &status_log = null_log()) const;

    const ECIContainer &eci() const;


  private:


    /// Pointer to PrimClex
    const PrimClex *m_primclex;

    /// Cluster expansion description: name, property, calctype, ref, bset, eci
    ClexDescription m_desc;

    /// Pointer to ClexBasis held by PrimClex
    mutable const ClexBasis *m_clex_basis;

    /// Clexulator (obtained from PrimClex, but must be separate object because
    /// of how Clexulator are evaluated)
    mutable Clexulator m_clexulator;

    /// Pointer to ECIContainer held by PrimClex
    mutable const ECIContainer *m_eci;

  };

  /// \brief Compare using descriptions: A.desc() < B.desc()
  bool operator<(const Clex &A, const Clex &B);

}

#endif
