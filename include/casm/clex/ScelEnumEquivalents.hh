#ifndef CASM_ScelEnumEquivalents
#define CASM_ScelEnumEquivalents

#include "casm/container/Array.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/symmetry/EnumEquivalents.hh"
#include "casm/clex/Supercell.hh"

namespace CASM {

  ENUMERATOR_TRAITS(ScelEnumEquivalents)

  class ScelEnumEquivalents :
    public EnumEquivalents <
    Supercell,
    Array<SymOp>::const_iterator,
    SymOp,
    SymRepIndexCompare > {

    // -- Required members -------------------

  public:

    ScelEnumEquivalents(const Supercell &scel);

    ENUMERATOR_MEMBERS(ScelEnumEquivalents)
  };

}

#endif
