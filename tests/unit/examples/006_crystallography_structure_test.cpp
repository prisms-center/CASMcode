#include "casm/crystallography/Structure.hh"
#include "crystallography/TestStructures.hh"  // for test::ZrO_prim
#include "gtest/gtest.h"

// Note: This example is still a stub
//
// A Structure is an object that holds a (shared, const) BasicStructure and its
// symmetry information. A Structure has:
// - a BasicStructure (std::shared_ptr<const xtal::BasicStructure>),
// representing the reference
//   infinite crystal and the DoF space
// - the basic structure's factor group (MasterSymGroup), which stores the
// symmetry operations
//   that leave the infinite crystal lattice, basis, and the DoF space
//   invariant, excluding operations that differ only by translational symmetry
// - keys (SymGroupRepID) used to access symmetry group representations
// (SymGroupRep) that are
//   stored in the MasterSymGroup object (accessed by
//   MasterSymGroup::representation). The SymGroupRep encode how symmetry
//   operations in the factor group do the following:
//   - permute sublattices ("basis_permutation_symrep_ID", SymGroupRepID)
//   - permute discrete site occupation values on each sublattice
//   ("occupant_symrep",
//     std::vector<SymGroupRepID>)
//   - transform global continuous DoF values ("global_dof_symrep_ID",
//   SymGroupRepID)
//   - transform local continuous DoF values on each sublattice
//   ("site_dof_symrep_IDs",
//     std::vector<std::map<DoFKey, SymGroupRepID>>)
//
// Note: In CASM, a Structure is often referred to as the "primitive crystal
// structure" or "prim", especially in the context where it is generally
// expected to the primitive unit cell.

TEST(ExampleCrystallographyStructure, StructureConstructor) {
  // A Structure is constructed from a BasicStructure.
  // - The BasicStructure held in Structure is const
  // - The factor group is constructed upon request
  //
  // In most contexts it is stored in a shared_ptr so that the BasicStructure
  // and factor group can
  //   be shared and used in various contexts.

  CASM::xtal::BasicStructure ZrO_basic_structure = test::ZrO_prim();
  auto ZrO_shared_prim =
      std::make_shared<const CASM::Structure>(ZrO_basic_structure);

  // TODO: more Structure examples
  EXPECT_EQ(ZrO_shared_prim->factor_group().size(), 24);
}
