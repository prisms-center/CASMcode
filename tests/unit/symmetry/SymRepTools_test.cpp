// This test is for SymRepTools, it would be cleanest if the setup and test
// cases were independent of Structure, clex, and enumerator, but they are
// currently used for setting up the SymGroupReps with known problematic cases

#include "casm/symmetry/io/json/SymRepTools.hh"

#include "autotools.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/DoFSet.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/enumerator/ConfigEnumInput_impl.hh"
#include "casm/enumerator/DoFSpace.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/symmetry/IrrepDecomposition.hh"
#include "casm/symmetry/IrrepDecompositionImpl.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymInfo.hh"
#include "casm/symmetry/io/json/SymGroup_json_io.hh"
#include "casm/symmetry/io/stream/SymInfo_stream_io.hh"
#include "gtest/gtest.h"

using namespace CASM;

namespace {
Eigen::Matrix3l _fcc_conventional_transf_mat() {
  Eigen::Matrix3l transf_mat;
  transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
  return transf_mat;
}

template <typename IrrepInfoVectorType>
void print_irreps(IrrepInfoVectorType const &irreps) {
  using SymRepTools_v2::IrrepDecompositionImpl::pretty;
  using SymRepTools_v2::IrrepDecompositionImpl::prettyc;

  std::cout << "irreps.size(): " << irreps.size() << std::endl;

  for (auto const &irrep : irreps) {
    std::cout << "---" << std::endl;
    std::cout << "symmetrized irrep: " << std::endl;
    std::cout << "index: " << irrep.index << std::endl;
    std::cout << "characters: " << prettyc(irrep.characters.transpose())
              << std::endl;
    std::cout << "subspace: \n" << prettyc(irrep.trans_mat) << std::endl;
    std::cout << "directions.size() (number of orbits): "
              << irrep.directions.size() << std::endl;
    Index orbit_index = 0;
    for (auto const &orbit : irrep.directions) {
      std::cout << "-" << std::endl;
      std::cout << "orbit: " << orbit_index << std::endl;
      for (auto const &direction : orbit) {
        std::cout << prettyc(direction.transpose()) << std::endl;
      }
      ++orbit_index;
    }
  }
  return;
}
}  // namespace

class IrrepDecompositionTest0 : public testing::Test {
 protected:
  // must make in SetUp
  std::shared_ptr<CASM::Structure const> shared_prim;
  std::shared_ptr<CASM::Supercell> shared_supercell;
  std::unique_ptr<DoFSpace> unique_dof_space;
  MasterSymGroup symrep_master_group;
  SymGroupRepID symrep_id;
  SymGroupRep const *symrep_ptr;

  IrrepDecompositionTest0() {}

  // Builds the collective DoF SymGroupRep which we want to test finding the
  // irreducible representations of
  void SetUp() override {
    // FCC, binary, with 3d displacements
    using namespace xtal;

    Molecule A = Molecule::make_atom("A");
    Molecule B = Molecule::make_atom("B");

    Lattice lat{Eigen::Vector3d{0.0, 2.0, 2.0}, Eigen::Vector3d{2.0, 0.0, 2.0},
                Eigen::Vector3d{2.0, 2.0, 0.0}};

    BasicStructure struc{lat};
    struc.set_basis({Site{Coordinate{0.0, 0.0, 0.0, lat, FRAC}, {A, B}}});
    struc.set_global_dofs({AnisoValTraits::strain("GL")});

    shared_prim = std::make_shared<CASM::Structure const>(struc);
    shared_supercell = std::make_shared<CASM::Supercell>(
        shared_prim, Eigen::Matrix3l::Identity());

    // Construct the disp DoF space.
    ConfigEnumInput config_input{*shared_supercell};
    DoFKey dof_key = "GLstrain";
    unique_dof_space =
        notstd::make_unique<DoFSpace>(make_dof_space(dof_key, config_input));
    SupercellSymInfo const &sym_info = shared_supercell->sym_info();
    std::vector<PermuteIterator> group = make_invariant_subgroup(config_input);
    EXPECT_EQ(group.size(), 48);

    symrep_ptr = &make_dof_space_symrep(*unique_dof_space, sym_info, group,
                                        symrep_master_group, symrep_id);
    EXPECT_EQ(symrep_ptr->dim(), 6);
  }

  void print_subspace() {
    Eigen::MatrixXd const &subspace = unique_dof_space->basis();
    std::cout << "subspace: " << std::endl;
    std::cout << subspace << std::endl << std::endl;
  }

  void print_symrep() {
    SymGroupRep const &rep = *symrep_ptr;
    std::cout << "collective_dof_symrep: " << std::endl;
    Index i = 0;
    for (CASM::SymOp const &op : symrep_master_group) {
      Eigen::MatrixXd matrix_rep = *(rep.MatrixXd(op));
      if (i == 0) {
        EXPECT_TRUE(matrix_rep.isIdentity(TOL));
      }
      std::cout << "---" << std::endl;
      std::cout << "i: " << i << "  op.index(): " << op.index() << std::endl;
      std::cout << "symop: "
                << brief_description(op, shared_prim->lattice(),
                                     SymInfoOptions{CART})
                << std::endl;
      std::cout << "matrix_rep: \n" << matrix_rep << std::endl;
      i++;
    }
  }
};

TEST_F(IrrepDecompositionTest0, Test1_a) {
  SymGroupRep const &rep = *symrep_ptr;
  Eigen::MatrixXd const &subspace = unique_dof_space->basis();

  bool calc_wedges = false;
  std::optional<VectorSpaceSymReport> symmetry_report;
  symmetry_report =
      vector_space_sym_report(rep, symrep_master_group, subspace, calc_wedges);
}

TEST_F(IrrepDecompositionTest0, Test1_b) {
  using namespace SymRepTools;

  using SymRepTools_v2::IrrepDecompositionImpl::pretty;
  using SymRepTools_v2::IrrepDecompositionImpl::prettyc;

  SymGroupRep const &rep = *symrep_ptr;
  SymGroup const &head_group = symrep_master_group;
  Eigen::MatrixXd const &subspace = unique_dof_space->basis();
  bool allow_complex = true;

  std::vector<SymRepTools::IrrepInfo> irreps =
      irrep_decomposition(rep, head_group, subspace, allow_complex);
  Eigen::MatrixXd symmetry_adapted_dof_subspace =
      full_trans_mat(irreps).adjoint();

  // print_irreps(irreps);
  // std::cout << "symmetry_adapted_subspace: \n"
  //           << pretty(symmetry_adapted_dof_subspace) << std::endl;
}

TEST_F(IrrepDecompositionTest0, Test2) {
  using namespace SymRepTools_v2;
  using SymRepTools_v2::IrrepDecompositionImpl::pretty;
  using SymRepTools_v2::IrrepDecompositionImpl::prettyc;

  SymGroupRep const &rep = *symrep_ptr;
  SymGroup const &head_group = symrep_master_group;
  Eigen::MatrixXd const &subspace = unique_dof_space->basis();

  MatrixRep matrix_rep;
  for (Index i = 0; i < rep.size(); ++i) {
    matrix_rep.push_back(*rep[i]->MatrixXd());
  }
  GroupIndices head_group_indices;
  for (SymOp const &op : head_group) {
    head_group_indices.insert(op.index());
  }
  GroupIndicesOrbitVector all_subgroups = head_group.subgroups();
  GroupIndicesOrbitVector cyclic_subgroups = head_group.small_subgroups();
  bool allow_complex = true;
  IrrepDecomposition irrep_decomposition{matrix_rep,    head_group_indices,
                                         subspace,      cyclic_subgroups,
                                         all_subgroups, allow_complex};

  // print_irreps(irrep_decomposition.irreps);
  // std::cout << "symmetry_adapted_subspace: \n"
  //           << pretty(irrep_decomposition.symmetry_adapted_subspace)
  //           << std::endl;
}

class IrrepDecompositionTest1 : public testing::Test {
 protected:
  // must make in SetUp
  std::shared_ptr<CASM::Structure const> shared_prim;
  std::shared_ptr<CASM::Supercell> shared_supercell;
  std::unique_ptr<DoFSpace> unique_dof_space;
  MasterSymGroup symrep_master_group;
  SymGroupRepID symrep_id;
  SymGroupRep const *symrep_ptr;

  IrrepDecompositionTest1() {}

  // Builds the collective DoF SymGroupRep which we want to test finding the
  // irreducible representations of
  void SetUp() override {
    // FCC, binary, with 3d displacements
    using namespace xtal;

    Molecule A = Molecule::make_atom("A");
    Molecule B = Molecule::make_atom("B");

    SiteDoFSet disp_xyz{
        AnisoValTraits::disp(),       // AnisoVal type
        {"d0x", "d0y", "d0z"},        // axes names
        Eigen::Matrix3d::Identity(),  // basis
        {}                            // excluded_occs
    };

    // fails:
    // Lattice lat{Eigen::Vector3d{2.0, 2.0, 0.0},
    // Eigen::Vector3d{0.0, 2.0, 2.0},
    //             Eigen::Vector3d{2.0, 0.0, 2.0}};
    // Eigen::Matrix3l T;
    // T << 0,  1, -2,
    //      0, 1, 2,
    //      1, -1, 0;

    // succeeds
    Lattice lat{Eigen::Vector3d{0.0, 2.0, 2.0}, Eigen::Vector3d{2.0, 0.0, 2.0},
                Eigen::Vector3d{2.0, 2.0, 0.0}};
    Eigen::Matrix3l T = _fcc_conventional_transf_mat();

    BasicStructure struc{lat};
    struc.set_basis(
        {Site{Coordinate{0.0, 0.0, 0.0, lat, FRAC}, {A, B}, {disp_xyz}}});

    shared_prim = std::make_shared<CASM::Structure const>(struc);
    shared_supercell = std::make_shared<CASM::Supercell>(shared_prim, T);

    // Construct the disp DoF space.
    ConfigEnumInput config_input{*shared_supercell};
    DoFKey dof_key = "disp";
    unique_dof_space =
        notstd::make_unique<DoFSpace>(make_dof_space(dof_key, config_input));
    SupercellSymInfo const &sym_info = shared_supercell->sym_info();
    std::vector<PermuteIterator> group = make_invariant_subgroup(config_input);
    EXPECT_EQ(group.size(), 48 * 4);

    symrep_ptr = &make_dof_space_symrep(*unique_dof_space, sym_info, group,
                                        symrep_master_group, symrep_id);
    EXPECT_EQ(symrep_ptr->dim(), 12);
  }

  void print_subspace() {
    Eigen::MatrixXd const &subspace = unique_dof_space->basis();
    std::cout << "subspace: " << std::endl;
    std::cout << subspace << std::endl << std::endl;
  }

  void print_symrep() {
    SymGroupRep const &rep = *symrep_ptr;
    std::cout << "collective_dof_symrep: " << std::endl;
    Index i = 0;
    for (CASM::SymOp const &op : symrep_master_group) {
      Eigen::MatrixXd matrix_rep = *(rep.MatrixXd(op));
      if (i == 0) {
        EXPECT_TRUE(matrix_rep.isIdentity(TOL));
      }
      std::cout << "---" << std::endl;
      std::cout << "i: " << i << "  op.index(): " << op.index() << std::endl;
      std::cout << "symop: "
                << brief_description(op, shared_prim->lattice(),
                                     SymInfoOptions{CART})
                << std::endl;
      std::cout << "matrix_rep: \n" << matrix_rep << std::endl;
      i++;
    }
  }
};

TEST_F(IrrepDecompositionTest1, Test1_a) {
  SymGroupRep const &rep = *symrep_ptr;
  Eigen::MatrixXd const &subspace = unique_dof_space->basis();

  bool calc_wedges = false;
  std::optional<VectorSpaceSymReport> symmetry_report;
  symmetry_report =
      vector_space_sym_report(rep, symrep_master_group, subspace, calc_wedges);
}

TEST_F(IrrepDecompositionTest1, Test1_b) {
  using namespace SymRepTools;
  using SymRepTools_v2::IrrepDecompositionImpl::pretty;
  using SymRepTools_v2::IrrepDecompositionImpl::prettyc;

  SymGroupRep const &rep = *symrep_ptr;
  SymGroup const &head_group = symrep_master_group;
  Eigen::MatrixXd const &subspace = unique_dof_space->basis();
  bool allow_complex = true;

  std::vector<SymRepTools::IrrepInfo> irreps =
      irrep_decomposition(rep, head_group, subspace, allow_complex);
  Eigen::MatrixXd symmetry_adapted_dof_subspace =
      full_trans_mat(irreps).adjoint();

  // print_irreps(irreps);
  // std::cout << "symmetry_adapted_subspace: \n"
  //           << pretty(symmetry_adapted_dof_subspace) << std::endl;
}

TEST_F(IrrepDecompositionTest1, Test2) {
  using namespace SymRepTools_v2;
  using SymRepTools_v2::IrrepDecompositionImpl::pretty;
  using SymRepTools_v2::IrrepDecompositionImpl::prettyc;

  SymGroupRep const &rep = *symrep_ptr;
  SymGroup const &head_group = symrep_master_group;
  Eigen::MatrixXd const &subspace = unique_dof_space->basis();

  MatrixRep matrix_rep;
  for (Index i = 0; i < rep.size(); ++i) {
    matrix_rep.push_back(*rep[i]->MatrixXd());
  }
  GroupIndices head_group_indices;
  for (SymOp const &op : head_group) {
    head_group_indices.insert(op.index());
  }
  GroupIndicesOrbitVector all_subgroups = head_group.subgroups();
  GroupIndicesOrbitVector cyclic_subgroups = head_group.small_subgroups();
  bool allow_complex = true;
  IrrepDecomposition irrep_decomposition{matrix_rep,    head_group_indices,
                                         subspace,      cyclic_subgroups,
                                         all_subgroups, allow_complex};

  // print_irreps(irrep_decomposition.irreps);
  // std::cout << "symmetry_adapted_subspace: \n"
  //           << pretty(irrep_decomposition.symmetry_adapted_subspace)
  //           << std::endl;
}

class IrrepDecompositionTest2 : public testing::Test {
 protected:
  // must make in SetUp
  std::shared_ptr<CASM::Structure const> shared_prim;
  std::shared_ptr<CASM::Supercell> shared_supercell;
  std::unique_ptr<DoFSpace> unique_dof_space;
  MasterSymGroup symrep_master_group;
  SymGroupRepID symrep_id;
  SymGroupRep const *symrep_ptr;

  IrrepDecompositionTest2() {}

  // Builds the collective DoF SymGroupRep which we want to test finding the
  // irreducible representations of
  void SetUp() override {
    // FCC base structure,
    // - with corner atoms {A, B} allowed to displace in 3d
    // - with face atoms {C, D} allowed to displace in the face plane,
    // -> test case where all atoms may displace, but along different dimensions
    //    so no rigid translations are possible
    using namespace xtal;

    Molecule A = Molecule::make_atom("A");
    Molecule B = Molecule::make_atom("B");
    Molecule C = Molecule::make_atom("C");
    Molecule D = Molecule::make_atom("D");

    SiteDoFSet disp_xyz{
        AnisoValTraits::disp(),       // AnisoVal type
        {"d0x", "d0y", "d0z"},        // axes names
        Eigen::Matrix3d::Identity(),  // basis
        {}                            // excluded_occs
    };

    SiteDoFSet disp_xy{
        AnisoValTraits::disp(),  // AnisoVal type
        {"d1x", "d1y"},          // axes names
        Eigen::MatrixXd({        // basis: allow displacements in xy
                         {1., 0., 0.},
                         {0., 1., 0.}})
            .transpose(),
        {}  // excluded_occs
    };

    SiteDoFSet disp_yz{
        AnisoValTraits::disp(),  // AnisoVal type
        {"d2y", "d2z"},          // axes names
        Eigen::MatrixXd({        // basis: allow displacements in yz
                         {0., 1., 0.},
                         {0., 0., 1.}})
            .transpose(),
        {}  // excluded_occs
    };

    SiteDoFSet disp_xz{
        AnisoValTraits::disp(),  // AnisoVal type
        {"d3x", "d3z"},          // axes names
        Eigen::MatrixXd({        // basis: allow displacements in xz
                         {1., 0., 0.},
                         {0., 0., 1.}})
            .transpose(),
        {}  // excluded_occs
    };

    Lattice lat{Eigen::Vector3d{4.0, 0.0, 0.0}, Eigen::Vector3d{0.0, 4.0, 0.0},
                Eigen::Vector3d{0.0, 0.0, 4.0}};

    BasicStructure struc{lat};
    struc.set_basis(
        {Site{Coordinate{0.0, 0.0, 0.0, lat, FRAC}, {A, B}, {disp_xyz}},
         Site{Coordinate{0.5, 0.5, 0.0, lat, FRAC}, {C, D}, {disp_xy}},
         Site{Coordinate{0.0, 0.5, 0.5, lat, FRAC}, {C, D}, {disp_yz}},
         Site{Coordinate{0.5, 0.0, 0.5, lat, FRAC}, {C, D}, {disp_xz}}});

    shared_prim = std::make_shared<CASM::Structure const>(struc);
    shared_supercell = std::make_shared<CASM::Supercell>(
        shared_prim, Eigen::Matrix3l::Identity());

    // std::cout << "prim factor group: " << std::endl;
    // SymInfoOptions opt{CART};
    // brief_description(log(), shared_prim->factor_group(),
    //                   shared_prim->lattice(), opt);

    // Construct the disp DoF space.
    ConfigEnumInput config_input{*shared_supercell};
    DoFKey dof_key = "disp";
    unique_dof_space =
        notstd::make_unique<DoFSpace>(make_dof_space(dof_key, config_input));
    SupercellSymInfo const &sym_info = shared_supercell->sym_info();
    std::vector<PermuteIterator> group = make_invariant_subgroup(config_input);
    EXPECT_EQ(group.size(), 48);

    symrep_ptr = &make_dof_space_symrep(*unique_dof_space, sym_info, group,
                                        symrep_master_group, symrep_id);
    EXPECT_EQ(symrep_ptr->dim(), 9);
    // print_symrep();
  }

  void print_subspace() {
    Eigen::MatrixXd const &subspace = unique_dof_space->basis();
    std::cout << "subspace: " << std::endl;
    std::cout << subspace << std::endl << std::endl;
  }

  void print_symrep() {
    SymGroupRep const &rep = *symrep_ptr;
    std::cout << "collective_dof_symrep: " << std::endl;
    Index i = 0;
    for (CASM::SymOp const &op : symrep_master_group) {
      Eigen::MatrixXd matrix_rep = *(rep.MatrixXd(op));
      if (i == 0) {
        EXPECT_TRUE(matrix_rep.isIdentity(TOL));
      }
      std::cout << "---" << std::endl;
      std::cout << "i: " << i << "  op.index(): " << op.index() << std::endl;
      std::cout << "symop: "
                << brief_description(op, shared_prim->lattice(),
                                     SymInfoOptions{CART})
                << std::endl;
      std::cout << "matrix_rep: \n" << matrix_rep << std::endl;
      i++;
    }
  }
};

TEST_F(IrrepDecompositionTest2, Test1_a) {
  SymGroupRep const &rep = *symrep_ptr;
  Eigen::MatrixXd const &subspace = unique_dof_space->basis();

  bool calc_wedges = false;
  std::optional<VectorSpaceSymReport> symmetry_report;
  symmetry_report =
      vector_space_sym_report(rep, symrep_master_group, subspace, calc_wedges);
}

TEST_F(IrrepDecompositionTest2, Test1_b) {
  using namespace SymRepTools;
  using SymRepTools_v2::IrrepDecompositionImpl::pretty;
  using SymRepTools_v2::IrrepDecompositionImpl::prettyc;

  SymGroupRep const &rep = *symrep_ptr;
  SymGroup const &head_group = symrep_master_group;
  Eigen::MatrixXd const &subspace = unique_dof_space->basis();
  bool allow_complex = true;

  std::vector<SymRepTools::IrrepInfo> irreps =
      irrep_decomposition(rep, head_group, subspace, allow_complex);
  Eigen::MatrixXd symmetry_adapted_dof_subspace =
      full_trans_mat(irreps).adjoint();

  // print_irreps(irreps);
  // std::cout << "symmetry_adapted_subspace: \n"
  //           << pretty(symmetry_adapted_dof_subspace) << std::endl;
}

TEST_F(IrrepDecompositionTest2, Test2) {
  using namespace SymRepTools_v2;
  using SymRepTools_v2::IrrepDecompositionImpl::pretty;
  using SymRepTools_v2::IrrepDecompositionImpl::prettyc;

  SymGroupRep const &rep = *symrep_ptr;
  SymGroup const &head_group = symrep_master_group;
  Eigen::MatrixXd const &subspace = unique_dof_space->basis();

  MatrixRep matrix_rep;
  for (Index i = 0; i < rep.size(); ++i) {
    matrix_rep.push_back(*rep[i]->MatrixXd());
  }
  GroupIndices head_group_indices;
  for (SymOp const &op : head_group) {
    head_group_indices.insert(op.index());
  }
  GroupIndicesOrbitVector all_subgroups = head_group.subgroups();
  GroupIndicesOrbitVector cyclic_subgroups = head_group.small_subgroups();
  bool allow_complex = true;
  IrrepDecomposition irrep_decomposition{matrix_rep,    head_group_indices,
                                         subspace,      cyclic_subgroups,
                                         all_subgroups, allow_complex};

  // print_irreps(irrep_decomposition.irreps);
  // std::cout << "symmetry_adapted_subspace: \n"
  //           << irrep_decomposition.symmetry_adapted_subspace << std::endl;
}
