#include <iterator>
#include "casm/crystallography/SymTools.hh"
#include "casm/symmetry/SymTools.hh"
#include "casm/symmetry/SymTools_impl.hh"
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
      auto subgroup_operation_indices = xtal::invariant_subgroup_indices(lat, super_group);
      return subgroup_from_indices(super_group, subgroup_operation_indices);
    }
  } // namespace sym
} // namespace CASM
