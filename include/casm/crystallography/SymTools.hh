#ifndef XTALSYMTOOLS_HH
#define XTALSYMTOOLS_HH

#include <vector>;
#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/LatticeIsEquivalent.hh"
#include "casm/CASM_global_definitions.hh"

namespace CASM {
  namespace xtal {
    class Lattice;

    /// \brief Construct indices of the subgroup that leaves a lattice unchanged
    std::vector<Index> invariant_subgroup_indices(const Lattice &lat, SymOpVector const &super_grp);

    template <typename ExternSymOpVector>
    std::vector<Index> invariant_subgroup_indices(const Lattice &lat, ExternSymOpVector const &super_grp) {
      return invariant_subgroup_indices(lat, adapter::Adapter<SymOpVector, ExternSymOpVector>()(super_grp));
    }

    /// \brief Construct indices of the subgroup for which this->is_equivalent(copy_apply(op, *this))
    template <typename OutputIt>
    OutputIt invariant_subgroup_indices(const Lattice &lat, const std::vector<SymOp> &super_group, OutputIt result) {
      LatticeIsEquivalent is_equiv(lat);

      Index ix = 0;
      for(auto it = super_group.begin(); it != super_group.end(); ++it) {
        if(is_equiv(*it)) {
          *result = ix;
          ++result;
        }
        ++ix;
      }
      return result;
    }
  } // namespace xtal
} // namespace CASM

#endif
