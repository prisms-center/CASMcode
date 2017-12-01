#include "casm/clex/ScelEnumEquivalents.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/Lattice_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/NeighborList.hh"

namespace CASM {

  namespace {

    struct MakeInvariantSubgroup {

      MakeInvariantSubgroup() {}

      template<typename SymOpIterator, typename SymOpOutputIterator>
      SymOpOutputIterator operator()(const Supercell &scel, SymOpIterator begin, SymOpIterator end, SymOpOutputIterator result) {
        return scel.lattice().invariant_subgroup(begin, end, result);
      }
    };

  }

  const std::string ScelEnumEquivalents::enumerator_name = "ScelEnumEquivalents";

  ScelEnumEquivalents::ScelEnumEquivalents(const Supercell &scel) :
    EnumEquivalents<Supercell, Array<SymOp>::const_iterator, SymOp, SymRepIndexCompare>(
      scel.canonical_form(),
      scel.prim().point_group().begin(),
      scel.prim().point_group().end(),
      MakeInvariantSubgroup()) {}

}

