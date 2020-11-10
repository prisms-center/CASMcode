#include "casm/clex/ScelEnumEquivalents.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/container/Permutation.hh"
#include "casm/symmetry/SymTools.hh"
#include "casm/symmetry/SymTools_impl.hh"

namespace CASM {

  namespace {

    struct MakeInvariantSubgroup {

      MakeInvariantSubgroup() {}

      template<typename SymOpIterator, typename SymOpOutputIterator>
      SymOpOutputIterator operator()(const Supercell &scel, SymOpIterator begin, SymOpIterator end, SymOpOutputIterator result) {
        return sym::invariant_subgroup(std::vector<SymOp>(begin, end), scel.lattice(), result);
      }
    };

  }

  const std::string ScelEnumEquivalents::enumerator_name = "ScelEnumEquivalents";

  ScelEnumEquivalents::ScelEnumEquivalents(const Supercell &scel) :
    EnumEquivalents<Supercell, std::vector<SymOp>::const_iterator, SymOp, SymRepIndexCompare>(
      scel.canonical_form(),
      scel.prim().point_group().begin(),
      scel.prim().point_group().end(),
      MakeInvariantSubgroup()) {}

}
