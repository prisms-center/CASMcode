#ifndef CASM_symmetry_LatticeEnumEquivalents
#define CASM_symmetry_LatticeEnumEquivalents

#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/symmetry/EnumEquivalents.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/symmetry/SymOpRepresentation.hh"

namespace CASM {
class SymGroup;

using sym::copy_apply;  // TODO: Base class for equivalents still relies on
                        // CASM::copy_apply.

/// \brief Enumerate equivalent Lattics, given a SymGroup
///
/// - The 'begin' Lattice is always the canonical form, with respect to the
///   specified SymGroup
/// - Currently requires super_g to have a valid MasterSymGroup. This
/// requirement could be relaxed
///   if the comparison function is changed.
///
/// \ingroup EnumEquivalents
/// \ingroup LatticeEnum
///
class LatticeEnumEquivalents
    : public EnumEquivalents<xtal::Lattice,
                             std::vector<CASM::SymOp>::const_iterator,
                             CASM::SymOp, SymRepIndexCompare> {
 public:
  LatticeEnumEquivalents(const xtal::Lattice &lat, const SymGroup &super_g);

  std::string name() const override { return enumerator_name; }

  static const std::string enumerator_name;
};
}  // namespace CASM

#endif
