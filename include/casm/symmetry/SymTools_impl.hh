#ifndef SYMTOOLS_IMPL_HH
#define SYMTOOLS_IMPL_HH

#include "casm/symmetry/SymGroup.hh"
#include "casm/crystallography/Lattice.hh"

namespace CASM {
  namespace sym {
    template <typename OutputIt>
    OutputIt invariant_subgroup(const std::vector<SymOp> &super_group,
                                const xtal::Lattice &lat,
                                OutputIt result) {
      auto subgroup_operation_indices = lat.invariant_subgroup_indices(super_group);
      for(auto ix : subgroup_operation_indices) {
        auto it = super_group.begin();
        std::advance(it, ix);
        *result = *it;
        ++result;
      }
      return result;
    }
  }
}

#endif
