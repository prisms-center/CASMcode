#include "casm/clex/ScelEnumEquivalents.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/NeighborList.hh"

namespace CASM {

  namespace {

    struct MakeInvariantSubgroup {

      MakeInvariantSubgroup() {}

      template<typename SymOpIterator, typename SymOpOutputIterator>
      SymOpOutputIterator operator()(const Supercell &scel, SymOpIterator begin, SymOpIterator end, SymOpOutputIterator result) {
        double tol = scel.primclex().crystallography_tol();
        return scel.lattice().find_invariant_subgroup(begin, end, result, tol);
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

