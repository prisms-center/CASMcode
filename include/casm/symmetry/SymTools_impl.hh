#ifndef SYMTOOLS_IMPL_HH
#define SYMTOOLS_IMPL_HH

#include "casm/symmetry/SymGroup.hh"
#include "casm/crystallography/Lattice.hh"

namespace CASM {
  namespace sym {
    template <typename OutputIt>
    OutputIt invariant_subgroup(std::vector<SymOp>::const_iterator begin,
                                std::vector<SymOp>::const_iterator end,
                                const xtal::Lattice &lat,
                                OutputIt result) {
      auto subgroup_operation_indices = lat.invariant_subgroup_indices(begin, end);
      for(auto ix : subgroup_operation_indices) {
        auto it = begin;
        std::advance(it, ix);
        *result = *it;
        ++result;
      }
      return result;
    }
  }
}

#endif
