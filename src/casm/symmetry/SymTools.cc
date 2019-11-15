#include <iterator>
#include "casm/crystallography/SymTools.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/SymTools.hh"
#include "casm/symmetry/SymTools_impl.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymBasisPermute.hh"
#include "casm/misc/CASM_Eigen_math.hh"

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

    //TODO: Do we keep passing by reference or do we want to change our ways here
    //and start passing by pointer?
    //it could just be:
    //apply_symmetry(op, my_ucc)  // returns new value
    //apply_symmetry(op, &my_ucc) // modifies value
    xtal::UnitCellCoord &apply(const CASM::SymOp &op, xtal::UnitCellCoord &mutating_ucc, const xtal::Structure &prim) {

      // transform using stored SymBasisPermute representation
      const SymBasisPermute &rep = *op.get_basis_permute_rep(prim.basis_permutation_symrep_ID());
      mutating_ucc._unitcell() = rep.matrix() * mutating_ucc.unitcell() + rep[mutating_ucc.sublattice()].unitcell();
      mutating_ucc._sublattice() = rep[mutating_ucc.sublattice()].sublattice();

      // additional translations (such as needed for supercell factor groups),
      // are stored in SymOp::integral_tau() (in cartesian coordinates)
      // this converts that to fractional coordinates and adds it to this->unitcell()
      mutating_ucc._unitcell() += lround(prim.lattice().inv_lat_column_mat() * op.integral_tau());

      return mutating_ucc;
    }

    xtal::UnitCellCoord copy_apply(const CASM::SymOp &op, const xtal::UnitCellCoord &reference_ucc, const xtal::Structure &prim)  {
      UnitCellCoord result(reference_ucc);
      sym::apply(op, result, prim);
      return result;
    }

  } // namespace sym
} // namespace CASM
