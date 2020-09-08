#include "gtest/gtest.h"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Structure.hh"
#include "crystallography/TestStructures.hh" // for test::ZrO_prim

// This example introduces Configuration and related classes

// Configuration
// -------------
//
// A Configuration represents a particular periodic perturbation of the infinite crystal within
//   the space of allowed perturbations defined by the BasicStructure. Configuration has:
// - a Supercell, representing the translational periodicity of the perturbation
// - a ConfigDoF, representing the values of DoF (local discrete and local continuous)
//   "within" the supercell (i.e. the translationally unique perturbation), plus the values of the
//   global DoF.
// - calculated properties (std::map<std::string, MappedProperties>), the values of properties
//   (i.e. energy, relaxation displacements, relaxation strain, etc.) that are dependent on the
//   values of the DoF. Properties are stored in a map with key "calctype" to allow for different
//   values of the properties depending on the calculation method.
//
// Supercell
// ---------
//
//   A Supercell provides access to:
//   - SupercellSymInfo, which provides symmetry representations which transform site indices,
//     and values of global and site DoF given the invariance of the supercell lattice vectors.
//   - a SuperNeighborList, which specifies, for each unit cell in the supercell, the indices of the
//     translationally equivalent neighboring sites in a local region surronding the reference cell.
//     This provides fast access of neighboring site DoF values when evaluating cluster expansion
//     basis functions for all Configuration with the same Supercell. SuperNeighborList is
//     described more completely in example 401_clex_neighbor_list_test.cpp.
//
//     SupercellSymInfo has:
//     - a xtal::Superlattice,
//     - the supercell factor group (SymGroup), the subgroup of the prim factor group that
//       leaves the super lattice vectors invariant.
//     - the translation permutations ("translation_permutations", std::vector<Permutation>), which
//       describe the unique ways primitive structure lattice vector translations permute DoF values
//       in a ConfigDoF associated with the Supercell. The number of translation permutations is
//       equal the to volume of the supercell w.r.t. the primitive structure unit cell.
//     - keys (SymGroupRepID) used to access symmetry group representations (SymGroupRep) that
//       encode how discrete DoF values and local and global continuous DoF values transform due
//       to supercell factor group operations.
//     - converters between UnitCellCood site indices and the corresponding linear index into the
//       ConfigDoF occupation value vector, or column index in the continuous site DoF value matrix
//
//
// ConfigDoF
// ---------
//
// ConfigDoF has:
// - the values of the discrete site DoF within the supercell ("occupation", Eigen::VectorXi)
//
//   Example: Occupation values, accessed via `Eigen::VectorXi const &ConfigDoF::occupation()`,
//   with integer value corresponding to which Molecule in the the a Site::occupant_dof vector
//   is occupying a particular site:
//
//       [<- sublattice 0 "occ" values -> | <- sublattice 1 "occ" values -> | ... ]
//
// - the values of the continuous site DoF within the supercell ("local_dofs",
//   std::map<DoFKey, LocalDoFContainerType>).
//
//   Example: Displacement values, with DoF basis equal to the standard basis (dx, dy, dz),
//   accessed via `Eigen::MatrixXd const &ConfigDoF::local_dofs("disp").values()`:
//
//       [<- sublattice 0 dx values -> | <- sublattice 1 dx values -> | ... ]
//       [<- sublattice 0 dy values -> | <- sublattice 1 dy values -> | ... ]
//       [<- sublattice 0 dz values -> | <- sublattice 1 dz values -> | ... ]
//
// - the values of the continuous global DoF ("global_dofs",
//   std::map<DoFKey, GlobalDoFContainerType>).
//
//   Example: GLstrain values, with DoF basis equal to the standard basis, accessed via
//   `Eigen::VectorXd const &ConfigDoF::global_dofs("GLstrain").values()`:
//
//       [e_xx, e_yy, e_zz, sqrt(2)*e_yz, sqrt(2)*e_xz, sqrt(2)*e_xy]
//
//
// MappedProperties
// ----------------
//
/// Associated with a particular Configuration as "calculated properties"
/// (std::map<std::string, MappedProperties>), MappedProperties are values of properties (i.e.
/// energy, relaxation displacements, relaxation strain, etc.) that are dependent on the values of
/// the DoF.
///
/// MappedProperties are "mapped" in the sense that they are not necessarily the raw output values
/// from a DFT calculation (though some properties are), but are determined by "mapping" the
/// final calculated structure's atomic positions and lattice vectors to the ideal configuration
/// they are most similar to. In the case of significant atomic and lattice relaxations this
/// allows associating calculated properties with the correct Configuration. If the relaxations
/// are too severe, there may be no Configuration that is an appropriate mapping. The user has
/// control of criteria used to judge whether any Configuration is an appropriate mapping and
/// which is the best mapping.
///
/// Note: Multiple instances of MappedProperties may be associated with a single Configuration.
/// This may occur if multiple calculation types are performed for a single Configuration. This
/// may also occur if multiple ideal Configuration are unstable and relax to the same other
/// Configuration.
///
/// MappedProperties has:
/// - "global" (std::map<std::string, Eigen::MatrixXd>), a map of property name to value for
///   global properties properties of the configuration such as energy, strain metrics, or lattice
///   vectors. Property names must follow CASM property naming conventions as documented for
///   AnisoValTraits. This map also holds any scalar properties such as energy or mapping cost
///   ("lattice_deformation_cost", "atomic_deformation_cost", "total_cost"), which can be queried
///   and set using member functions.
/// - "site" (std::map<std::string, Eigen::MatrixXd>), a map of property name to value for
///   properties of particular sites such as displacements, forces, or atomic coordinates.
///   Property names must follow CASM property naming conventions as documented for AnisoValTraits.
///
///

TEST(ExampleEnumerationSupercell, SupercellConstructor) {

  // Construct the "prim" Structure, representing the infinite crystal and allowed perturbations.
  auto shared_prim = std::make_shared<CASM::Structure const>(test::ZrO_prim());

  // The prim to super lattice transformation matrix, T:
  //   super_lattice_column_matrix = T * prim_lattice_column_matrix
  Eigen::Matrix3i T;
  T << 2, 0, 0,
  0, 2, 0,
  0, 0, 2;

  // Construct the Supercell
  CASM::Supercell supercell {shared_prim, T};

  EXPECT_EQ(supercell.volume(), 8);
}

TEST(ExampleEnumerationConfiguration, ConfigurationConstructor) {

  // Construct the (shared) "prim" Structure, representing the infinite crystal and allowed perturbations.
  auto shared_prim = std::make_shared<CASM::Structure const>(test::ZrO_prim());

  // The prim to super lattice transformation matrix, T:
  //   super_lattice_column_matrix = T * prim_lattice_column_matrix
  Eigen::Matrix3i T;
  T << 2, 0, 0,
  0, 2, 0,
  0, 0, 2;

  // Construct the (shared) Supercell
  auto shared_supercell = std::make_shared<CASM::Supercell>(shared_prim, T);

  // Construct a Configuration
  CASM::Configuration config {shared_supercell};

  EXPECT_EQ(shared_supercell->volume(), 8);
}
