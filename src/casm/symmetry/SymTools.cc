#include <iterator>
#include "casm/symmetry/SymTools.hh"
#include "casm/symmetry/SymTools_impl.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/symmetry/SymGroup.hh"
namespace {
  CASM::SymGroup subgroup_from_indices(const CASM::SymGroup &super_group, const std::vector<CASM::Index> &subgroup_indices) {
    std::vector<CASM::SymOp> subgroup_operations;
    for(auto ix : subgroup_indices) {
      subgroup_operations.push_back(super_group[ix]);
    }

    return CASM::SymGroup(subgroup_operations, &(super_group.lattice()), super_group.periodicity());
  }

  //*******************************************************************************************
}
namespace CASM {
  namespace sym {
    SymGroup invariant_subgroup(const SymGroup &super_group, const xtal::Lattice &lat) {
      auto subgroup_operation_indices = lat.invariant_subgroup_indices(super_group);
      return subgroup_from_indices(super_group, subgroup_operation_indices);
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
