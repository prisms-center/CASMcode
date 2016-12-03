#ifndef CASM_ScelEnumEquivalents
#define CASM_ScelEnumEquivalents

#include "casm/container/Array.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/symmetry/EnumEquivalents.hh"
#include "casm/clex/Supercell.hh"

namespace CASM {

  /// \brief Enumerate equivalent Supercell
  ///
  /// \ingroup EnumEquivalents
  /// \ingroup ScelEnumGroup
  ///
  class ScelEnumEquivalents :
    public EnumEquivalents <
    Supercell,
    Array<SymOp>::const_iterator,
    SymOp,
    SymRepIndexCompare > {

    // -- Required members -------------------

  public:

    ScelEnumEquivalents(const Supercell &scel);

    std::string name() const override {
      return enumerator_name;
    }

    static const std::string enumerator_name;

  };

}

#endif
