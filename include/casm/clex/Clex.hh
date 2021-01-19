#ifndef CASM_Clex
#define CASM_Clex

#include "casm/clex/Clexulator.hh"
#include "casm/clex/ECIContainer.hh"

namespace CASM {

/** \defgroup ClexClex Clex
 *  \ingroup Clex
 *  \brief Relates to cluster expansion basis sets and ECI
 *  @{
 */

/// Pair of Clexulator and ECIContainer
struct Clex {
  Clexulator clexulator;
  ECIContainer eci;

  Clex(Clexulator const &_clexulator, ECIContainer const &_eci)
      : clexulator(_clexulator), eci(_eci) {}
};

/** @} */
}  // namespace CASM

#endif
