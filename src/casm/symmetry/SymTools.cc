#include <iterator>
#include "casm/symmetry/SymTools.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/symmetry/SymGroup.hh"

namespace CASM {
  namespace sym {
    SymGroup invariant_subgroup(const SymGroup &super_group, const xtal::Lattice &lat) {
      auto subgroup_operation_indices = lat.invariant_subgroup_indices(super_group);
      return super_group.subgroup_from_indices(subgroup_operation_indices);
    }

    std::vector<SymOp> invariant_subgroup(std::vector<SymOp>::const_iterator begin,
                                          std::vector<SymOp>::const_iterator end,
                                          const xtal::Lattice &lat) {
      auto subgroup_operation_indices = lat.invariant_subgroup_indices(begin, end);
      std::vector<SymOp> result;
      invariant_subgroup(begin, end, lat, std::back_inserter(result));
      return result;

    }

  } // namespace sym
} // namespace CASM
