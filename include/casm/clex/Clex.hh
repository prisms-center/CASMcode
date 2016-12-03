#ifndef CASM_Clex
#define CASM_Clex
#include <cstddef>

#include "casm/CASM_global_definitions.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/app/ProjectSettings.hh"

namespace CASM {

  class PrimClex;
  class ECIContainer;
  template<typename ClustType> class GenericOrbitree;
  class SiteCluster;
  typedef GenericOrbitree<SiteCluster> SiteOrbitree;

  /** \defgroup ClexClex Clex
   *  \ingroup Clex
   *  \brief Relates to cluster expansion basis sets and ECI
   *  @{
   */

  /// \brief Data structure used for cluster expansions
  class Clex {

  public:

    Clex() {}

    Clex(const PrimClex &_primclex) :
      m_primclex(&_primclex), m_orbitree(nullptr), m_eci(nullptr) {}

    Clex(const PrimClex &_primclex, const ClexDescription &_desc) :
      m_primclex(&_primclex), m_desc(_desc), m_orbitree(nullptr), m_eci(nullptr) {}


    const PrimClex &primclex() const;

    const ClexDescription &desc() const;

    const SiteOrbitree &orbitree() const;

    Clexulator &clexulator(Log &status_log = null_log()) const;

    const ECIContainer &eci() const;


  private:


    /// Pointer to PrimClex
    const PrimClex *m_primclex;

    /// Cluster expansion description: name, property, calctype, ref, bset, eci
    ClexDescription m_desc;

    /// Pointer to SiteOrbitree held by PrimClex
    mutable const SiteOrbitree *m_orbitree;

    /// Clexulator (obtained from PrimClex, but must be separate object because
    /// of how Clexulator are evaluated)
    mutable Clexulator m_clexulator;

    /// Pointer to ECIContainer held by PrimClex
    mutable const ECIContainer *m_eci;

  };

  /// \brief Compare using descriptions: A.desc() < B.desc()
  bool operator<(const Clex &A, const Clex &B);

  /** @} */
}

#endif
