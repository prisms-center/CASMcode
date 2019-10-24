#include "casm/crystallography/SymTools.hh"
namespace CASM {
  namespace xtal {
    /// \brief Construct the subgroup that leaves a lattice unchanged
    std::vector<Index> invariant_subgroup_indices(const Lattice &lat, std::vector<SymOp> const &super_grp) {
      std::vector<Index> result;
      invariant_subgroup_indices(lat, super_grp, std::back_inserter(result));
      return result;
    }

  } // namespace xtal
} // namespace CASM
