#ifndef CASM_LatticeEnumEquivalents
#define CASM_LatticeEnumEquivalents

#include "casm/symmetry/EnumEquivalents.hh"
#include "casm/symmetry/SymOpRepresentation.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/SymTools.hh"

//TODO:
//Why does this class exist again? Literally nothing is using it.
//It depends fundamentally on CASM::SymOp, so it's getting the boot.
//Don't keep it in xtal.

namespace CASM {
  class SymGroup;

  using sym::copy_apply;  //TODO: Base class for equivalents still relies on CASM::copy_apply.

  /// \brief Enumerate equivalent Lattics, given a SymGroup
  ///
  /// - The 'begin' Lattice is always the canonical form, with respect to the
  ///   specified SymGroup
  /// - Currently requires super_g to have a valid MasterSymGroup. This requirement could be relaxed
  ///   if the comparison function is changed.
  ///
  /// \ingroup EnumEquivalents
  /// \ingroup LatticeEnum
  ///
  class LatticeEnumEquivalents :
    public EnumEquivalents<Lattice, std::vector<CASM::SymOp>::const_iterator, CASM::SymOp, SymRepIndexCompare> {

  public:
    LatticeEnumEquivalents(const Lattice &lat, const SymGroup &super_g);

    std::string name() const override {
      return enumerator_name;
    }

    static const std::string enumerator_name;
  };
}

#endif
