#ifndef CASM_LatticeEnumEquivalents
#define CASM_LatticeEnumEquivalents

#include "casm/symmetry/EnumEquivalents.hh"
#include "casm/symmetry/SymOpRepresentation.hh"
#include "casm/crystallography/Lattice.hh"

namespace CASM {

  class SymGroup;

  /// \brief Enumerate equivalent Lattics, given a SymGroup
  ///
  /// - The 'begin' Lattice is always the canonical form, with respect to the
  ///   specified SymGroup
  ///
  /// \ingroup EnumEquivalents
  /// \ingroup LatticeEnum
  ///
  class LatticeEnumEquivalents :
    public EnumEquivalents<Lattice, Array<SymOp>::const_iterator, SymOp, SymRepIndexCompare> {

  public:
    LatticeEnumEquivalents(const Lattice &lat, const SymGroup &super_g, double tol);

    std::string name() const override {
      return enumerator_name;
    }

    static const std::string enumerator_name;
  };

}

#endif
