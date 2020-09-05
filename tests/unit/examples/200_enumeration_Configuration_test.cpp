#include "gtest/gtest.h"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Structure.hh"
#include "crystallography/TestStructures.hh" // for test::ZrO_prim

// This example introduces the Configuration and Supercell classes

// A Configuration represents a particular periodic perturbation of the infinite crystal within
//   the space of allowed perturbations defined by the BasicStructure. Configuration has:
// - a ConfigDoF, representing the values of DoF (local discrete and local continuous)
//   within the supercell (i.e. the perturbation within the periodically tiled unit), plus the
//   values of the global DoF.
// - a Supercell, representing the translational periodicity of the perturbation
//
//   ConfigDoF has:
//   - the values of the discrete site DoF within the supercell ("occupation", Eigen::VectorXi)
//
//     Example: Occupation values, accessed via `Eigen::VectorXi const &ConfigDoF::occupation()`,
//     with integer value corresponding to which Molecule in the the a Site::occupant_dof vector
//     is occupying a particular site:
//
//         [<- sublattice 0 "occ" values -> | <- sublattice 1 "occ" values -> | ... ]
//
//   - the values of the continuous site DoF within the supercell ("local_dofs",
//     std::map<DoFKey, LocalDoFContainerType>).
//
//     Example: Displacement values, with DoF basis equal to the standard basis (dx, dy, dz),
//     accessed via `Eigen::MatrixXd const &ConfigDoF::local_dofs("disp").values()`:
//
//         [<- sublattice 0 dx values -> | <- sublattice 1 dx values -> | ... ]
//         [<- sublattice 0 dy values -> | <- sublattice 1 dy values -> | ... ]
//         [<- sublattice 0 dz values -> | <- sublattice 1 dz values -> | ... ]
//
//   - the values of the continuous global DoF ("global_dofs",
//     std::map<DoFKey, GlobalDoFContainerType>).
//
//     Example: GLstrain values, with DoF basis equal to the standard basis, accessed via
//     `Eigen::VectorXd const &ConfigDoF::global_dofs("GLstrain").values()`:
//
//         [e_xx, e_yy, e_zz, sqrt(2)*e_yz, sqrt(2)*e_xz, sqrt(2)*e_xy]
//
// A Supercell provides access to:
// - SupercellSymInfo, which provides symmetry representations which transform site indices,
//   and values of global and site DoF given the invariance of the supercell lattice vectors.
// - a SuperNeighborList, which specifies, for each unit cell in the supercell, the indices of the
//   translationally equivalent neighboring sites in a local region surronding the reference cell.
//   This provides fast access of neighboring site DoF values when evaluating cluster expansion
//   basis functions for all Configuration with the same Supercell.
//
//   SupercellSymInfo has:
//   - a xtal::Superlattice,
//   - the supercell factor group (SymGroup), the subgroup of the prim factor group that
//     leaves the super lattice vectors invariant.
//   - the translation permutations ("translation_permutations", std::vector<Permutation>), which
//     describe the unique ways primitive structure lattice vector translations permute DoF values
//     in a ConfigDoF associated with the Supercell. The number of translation permutations is
//     equal the to volume of the supercell w.r.t. the primitive structure unit cell.
//   - keys (SymGroupRepID) used to access symmetry group representations (SymGroupRep) that
//     encode how discrete DoF values and local and global continuous DoF values transform due
//     to supercell factor group operations.
//   - converters between UnitCellCood site indices and the corresponding linear index into the
//     ConfigDoF occupation value vector, or column index in the continuous site DoF value matrix

TEST(SupercellTest, SupercellConstructor) {

  // Construct the "prim" Structure, representing the infinite crystal and allowed perturbations.
  auto shared_prim = std::make_shared<CASM::Structure const>(test::ZrO_prim());

  // The prim to super lattice transformation matrix, T:
  //   super_lattice_column_matrix = T * prim_lattice_column_matrix
  Eigen::Matrix3i T;
  T << 2, 0, 0,
  0, 2, 0,
  0, 0, 2;

  // Construct the Supercell
  Supercell supercell {shared_prim, T};

  EXPECT_EQ(supercell.volume(), 8);
}
