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

    template <typename OutputIt>
    OutputIt invariant_subgroup(std::vector<SymOp>::const_iterator begin,
                                std::vector<SymOp>::const_iterator end,
                                const xtal::Lattice &lat,
                                OutputIt result) {
      auto subgroup_operation_indices = lat.invariant_subgroup_indices(begin, end);
      for(auto ix : subgroup_operation_indices) {
        auto it = begin;
        std::advance(begin, ix);
        *result = *it;
        ++result;
      }
      return result;
    }
  } // namespace sym
} // namespace CASM
