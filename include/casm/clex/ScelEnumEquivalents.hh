#ifndef CASM_ScelEnumEquivalents
#define CASM_ScelEnumEquivalents

#include "casm/clex/Supercell.hh"
#include "casm/container/Array.hh"
#include "casm/symmetry/EnumEquivalents.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {

/// \brief Enumerate equivalent Supercell
///
/// \ingroup EnumEquivalents
/// \ingroup ScelEnumGroup
///
class ScelEnumEquivalents
    : public EnumEquivalents<Supercell, std::vector<SymOp>::const_iterator,
                             SymOp, SymRepIndexCompare> {
  // -- Required members -------------------

 public:
  ScelEnumEquivalents(const Supercell &scel);

  std::string name() const override { return enumerator_name; }

  static const std::string enumerator_name;
};

}  // namespace CASM

#endif
