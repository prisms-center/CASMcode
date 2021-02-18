#include "casm/enumerator/DoFSpace.hh"

#include "autotools.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/DoFSet.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/enumerator/ConfigEnumInput_impl.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymInfo.hh"
#include "casm/symmetry/io/json/SymGroup_json_io.hh"
#include "casm/symmetry/io/json/SymRepTools.hh"
#include "casm/symmetry/io/stream/SymInfo_stream_io.hh"
#include "crystallography/TestStructures.hh"
#include "gtest/gtest.h"

using namespace CASM;
using namespace test;

namespace {

Eigen::Matrix3l _fcc_conventional_transf_mat() {
  Eigen::Matrix3l transf_mat;
  transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
  return transf_mat;
}

// print local DoF symrep matrices to log
void print_local_dof_symreps(SupercellSymInfo const &sym_info) {
  jsonParser rep_json;
  for (auto const &pair : sym_info.local_dof_symreps()) {
    rep_json[pair.first] =
        jsonParser::array(pair.second.size(), jsonParser::object());
    Index b = 0;
    for (auto const &group_handle : pair.second) {
      write_matrix_rep(group_handle, rep_json[pair.first][b]);
      ++b;
    }
  }
  log() << rep_json << std::endl;
}

}  // namespace

class DoFSpaceTest : public testing::Test {
 protected:
  std::shared_ptr<CASM::Structure const> shared_prim;
  std::shared_ptr<CASM::Supercell> shared_supercell;  // conventional unit cell

  DoFSpaceTest()
      : shared_prim(std::make_shared<CASM::Structure const>(
            test::FCC_ternary_GLstrain_disp_prim())),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, _fcc_conventional_transf_mat())) {}
};

TEST_F(DoFSpaceTest, ConstructorTest1) {
  // Construct the GLstrain DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "GLstrain";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);
  EXPECT_EQ(dof_space.basis().rows(), 6);
}

TEST_F(DoFSpaceTest, ConstructorTest2) {
  // Construct the disp DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);
  EXPECT_EQ(dof_space.basis().rows(), 4 * 3);
}

TEST_F(DoFSpaceTest, JsonIOTest1) {
  // Construct the GLstrain DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "GLstrain";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);

  jsonParser dof_space_json;
  to_json(dof_space, dof_space_json, shared_supercell->name() + "/0");

  // // Uncomment to print dof_space:
  // log() << dof_space_json << std::endl;

  // a couple tests to check JSON was constructed successfully, but does not
  // check each value
  EXPECT_EQ(dof_space_json.size(), 8);
  EXPECT_EQ(dof_space_json.contains("dof"), true);
  EXPECT_EQ(dof_space_json["transformation_matrix_to_supercell"].is_null(),
            true);
  EXPECT_EQ(dof_space_json["sites"].is_null(), true);
  EXPECT_EQ(dof_space_json.contains("basis"), true);
  EXPECT_EQ(dof_space_json.contains("glossary"), true);
  EXPECT_EQ(dof_space_json["axis_site_index"].is_null(), true);
  EXPECT_EQ(dof_space_json["axis_dof_component"].is_null(), true);
  EXPECT_EQ(dof_space_json.contains("identifier"), true);
  EXPECT_EQ(dof_space_json.contains("state"), false);
}

TEST_F(DoFSpaceTest, JsonIOTest2) {
  // Construct the disp DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);

  jsonParser dof_space_json;
  to_json(dof_space, dof_space_json, shared_supercell->name() + "/0");

  // Uncomment to print dof_space:
  log() << dof_space_json << std::endl;

  // a couple tests to check JSON was constructed successfully, but does not
  // check each value
  EXPECT_EQ(dof_space_json.size(), 8);
  EXPECT_EQ(dof_space_json.contains("dof"), true);
  EXPECT_EQ(dof_space_json["transformation_matrix_to_supercell"].is_null(),
            false);
  EXPECT_EQ(dof_space_json["sites"].is_null(), false);
  EXPECT_EQ(dof_space_json.contains("basis"), true);
  EXPECT_EQ(dof_space_json.contains("glossary"), true);
  EXPECT_EQ(dof_space_json["axis_site_index"].is_null(), false);
  EXPECT_EQ(dof_space_json["axis_dof_component"].is_null(), false);
  EXPECT_EQ(dof_space_json.contains("identifier"), true);
  EXPECT_EQ(dof_space_json.contains("state"), false);
}

TEST_F(DoFSpaceTest, VectorSpaceSymReportTest1) {
  // Construct the GLstrain DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "GLstrain";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);

  // Construct the VectorSpaceSymReport (calc_wedges==false)
  auto const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  bool calc_wedges = false;
  VectorSpaceSymReport report = vector_space_sym_report(
      dof_space, sym_info, invariant_group, calc_wedges);

  jsonParser report_json;
  to_json(report, report_json);

  // // Uncomment to print VectorSpaceSymReport:
  // log() << report_json << std::endl;

  EXPECT_EQ(report.symgroup_rep.size(), 48);
  EXPECT_EQ(report.irreps.size(), 3);
  EXPECT_EQ(report.irreducible_wedge.size(), 0);
  EXPECT_EQ(report.symmetry_adapted_dof_subspace.rows(), 6);
  EXPECT_EQ(report.symmetry_adapted_dof_subspace.cols(), 6);
  EXPECT_EQ(report.axis_glossary.size(), 6);
}

TEST_F(DoFSpaceTest, VectorSpaceSymReportTest2) {
  // Construct the GLstrain DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "GLstrain";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);

  // Construct the VectorSpaceSymReport (calc_wedges==true)
  auto const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  bool calc_wedges = true;
  VectorSpaceSymReport report = vector_space_sym_report(
      dof_space, sym_info, invariant_group, calc_wedges);
  jsonParser report_json;
  to_json(report, report_json);

  // // Uncomment to print VectorSpaceSymReport:
  // log() << report_json << std::endl;

  EXPECT_EQ(report.symgroup_rep.size(), 48);
  EXPECT_EQ(report.irreps.size(), 3);
  EXPECT_EQ(report.irreducible_wedge.size(), 6);
  EXPECT_EQ(report.symmetry_adapted_dof_subspace.rows(), 6);
  EXPECT_EQ(report.symmetry_adapted_dof_subspace.cols(), 6);
  EXPECT_EQ(report.axis_glossary.size(), 6);
}

TEST_F(DoFSpaceTest, VectorSpaceSymReportTest3) {
  // Construct the GLstrain DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);

  // Construct the VectorSpaceSymReport (calc_wedges==false)
  auto const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  bool calc_wedges = false;
  VectorSpaceSymReport report = vector_space_sym_report(
      dof_space, sym_info, invariant_group, calc_wedges);

  jsonParser report_json;
  to_json(report, report_json);

  // // Uncomment to print VectorSpaceSymReport:
  // log() << report_json << std::endl;

  EXPECT_EQ(report.symgroup_rep.size(), 192);
  EXPECT_EQ(report.irreps.size(), 3);
  EXPECT_EQ(report.irreducible_wedge.size(), 0);
  EXPECT_EQ(report.symmetry_adapted_dof_subspace.rows(), 12);
  EXPECT_EQ(report.symmetry_adapted_dof_subspace.cols(), 12);
  EXPECT_EQ(report.axis_glossary.size(), 12);
}

TEST_F(DoFSpaceTest, VectorSpaceSymReportTest4) {
  // Construct the GLstrain DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);

  // Construct the VectorSpaceSymReport (calc_wedges==true)
  auto const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  bool calc_wedges = true;
  VectorSpaceSymReport report = vector_space_sym_report(
      dof_space, sym_info, invariant_group, calc_wedges);
  jsonParser report_json;
  to_json(report, report_json);

  // // Uncomment to print VectorSpaceSymReport:
  // log() << report_json << std::endl;

  EXPECT_EQ(report.symgroup_rep.size(), 192);
  EXPECT_EQ(report.irreps.size(), 3);
  EXPECT_EQ(report.irreducible_wedge.size(), 2304);
  EXPECT_EQ(report.symmetry_adapted_dof_subspace.rows(), 12);
  EXPECT_EQ(report.symmetry_adapted_dof_subspace.cols(), 12);
  EXPECT_EQ(report.axis_glossary.size(), 12);
}

TEST_F(DoFSpaceTest, ExcludeHomogeneousModeSpace) {
  // Construct the disp DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space_0 = make_dof_space(dof_key, config_input);

  Eigen::MatrixXd homogeneous_mode_space =
      make_homogeneous_mode_space(dof_space_0);
  EXPECT_EQ(homogeneous_mode_space.rows(), dof_space_0.basis().rows());
  EXPECT_EQ(homogeneous_mode_space.cols(), 3);

  DoFSpace dof_space_1 = exclude_homogeneous_mode_space(dof_space_0);
  EXPECT_EQ(dof_space_1.basis().rows(), dof_space_0.basis().rows());
  EXPECT_EQ(dof_space_1.basis().cols(), dof_space_0.basis().cols() - 3);

  SupercellSymInfo const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  bool calc_wedges = false;
  std::optional<VectorSpaceSymReport> sym_report;
  DoFSpace dof_space_2 = make_symmetry_adapted_dof_space(
      dof_space_1, sym_info, invariant_group, calc_wedges, sym_report);
  EXPECT_EQ(dof_space_2.basis().rows(), dof_space_1.basis().rows());
  EXPECT_EQ(dof_space_2.basis().cols(), dof_space_1.basis().cols());

  // // check symmetry report
  // jsonParser json;
  // to_json(dof_space_2, json, "test", config_input, sym_report);
  // std::cout << json << std::endl;
}

class RestrictedLocalDoFSpaceTest : public testing::Test {
 protected:
  std::shared_ptr<CASM::Structure const> shared_prim;
  std::shared_ptr<CASM::Supercell> shared_supercell;  // conventional unit cell

  static xtal::BasicStructure make_prim();

  RestrictedLocalDoFSpaceTest()
      : shared_prim(std::make_shared<CASM::Structure const>(make_prim())),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, _fcc_conventional_transf_mat())) {}
};

xtal::BasicStructure RestrictedLocalDoFSpaceTest::make_prim() {
  using namespace xtal;

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  SiteDoFSet disp_xy{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d1", "d2"},            // axes names
      Eigen::MatrixXd({        // basis: allow displacements in xy
                       {1., 0., 0.},
                       {0., 1., 0.}})
          .transpose(),
      {}  // excluded_occs
  };

  Lattice lat{Eigen::Vector3d{0.000000000000, 1.754750223661, 1.754750223661},
              Eigen::Vector3d{1.754750223661, 0.000000000000, 1.754750223661},
              Eigen::Vector3d{1.754750223661, 1.754750223661, 0.000000000000}};

  BasicStructure struc{lat};
  struc.set_basis(
      {Site{Coordinate{0.0, 0.0, 0.0, lat, FRAC}, {A, B}, {disp_xy}}});
  return struc;
}

TEST_F(RestrictedLocalDoFSpaceTest, FactorGroupSize) {
  auto const &factor_group = shared_prim->factor_group();

  SymInfoOptions opt{CART};
  brief_description(log(), factor_group, shared_prim->lattice(), opt);

  EXPECT_EQ(factor_group.size(), 16);
}

TEST_F(RestrictedLocalDoFSpaceTest, PrimSiteDoFSymReps) {
  auto const &factor_group = shared_prim->factor_group();

  // check dimensions of prim disp symrep matrices, should be 2x2
  for (auto const &site_symrep_IDs : shared_prim->site_dof_symrep_IDs()) {
    SymGroupRep const &rep =
        factor_group.representation(site_symrep_IDs.find("disp")->second);
    for (SymOp const &op : factor_group) {
      Eigen::MatrixXd symrep_matrix = *(rep.MatrixXd(op));
      EXPECT_EQ(symrep_matrix.rows(), 2);
      EXPECT_EQ(symrep_matrix.cols(), 2);
    }
  }
}

TEST_F(RestrictedLocalDoFSpaceTest, CollectiveDoFSymReps) {
  // check dimensions of collective dof symrep matrices:
  // for all sites, should be 8x8

  // Construct the disp DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);
  SupercellSymInfo const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  EXPECT_EQ(invariant_group.size(), 64);

  MasterSymGroup symrep_master_group;
  SymGroupRepID id;
  SymGroupRep const &rep = make_dof_space_symrep(
      dof_space, sym_info, invariant_group, symrep_master_group, id);
  EXPECT_EQ(rep.size(), 64);

  std::cout << "collective_dof_symrep" << std::endl;
  Index i = 0;
  for (SymOp const &op : symrep_master_group) {
    Eigen::MatrixXd matrix_rep = *(rep.MatrixXd(op));
    if (i == 0) {
      EXPECT_TRUE(matrix_rep.isIdentity(TOL));
    }
    std::cout << "---" << std::endl;
    std::cout << "i: " << i << "  op.index(): " << op.index() << std::endl;
    std::cout << matrix_rep << std::endl;
    EXPECT_EQ(matrix_rep.rows(), 8);
    EXPECT_EQ(matrix_rep.cols(), 8);
    i++;
  }
}

TEST_F(RestrictedLocalDoFSpaceTest, SubsetCollectiveDoFSymReps) {
  // check dimensions of collective dof symrep matrices:
  // for sites {2, 3}, should be 4x4

  // Construct the disp DoF space, but only on sites 2 and 3
  ConfigEnumInput config_input{*shared_supercell, {2, 3}};
  DoFKey dof_key = "disp";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);
  SupercellSymInfo const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  EXPECT_EQ(invariant_group.size(), 16);

  MasterSymGroup symrep_master_group;
  SymGroupRepID id;
  SymGroupRep const &rep = make_dof_space_symrep(
      dof_space, sym_info, invariant_group, symrep_master_group, id);
  EXPECT_EQ(rep.size(), 16);

  std::cout << "collective_dof_symrep" << std::endl;
  Index i = 0;
  for (SymOp const &op : symrep_master_group) {
    Eigen::MatrixXd matrix_rep = *(rep.MatrixXd(op));
    if (i == 0) {
      EXPECT_TRUE(matrix_rep.isIdentity(TOL));
    }
    std::cout << "---" << std::endl;
    std::cout << "i: " << i << "  op.index(): " << op.index() << std::endl;
    std::cout << matrix_rep << std::endl;
    EXPECT_EQ(matrix_rep.rows(), 4);
    EXPECT_EQ(matrix_rep.cols(), 4);
    i++;
  }
}

TEST_F(RestrictedLocalDoFSpaceTest, ExcludeHomogeneousModeSpace) {
  // print_local_dof_symreps(shared_supercell->sym_info());

  // Construct the restricted disp DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space_0 = make_dof_space(dof_key, config_input);
  // std::cout << "including homogeneous_mode_space: \n"
  //           << dof_space_0.basis() << std::endl;
  EXPECT_EQ(dof_space_0.basis().rows(), 8);
  EXPECT_EQ(dof_space_0.basis().cols(), 8);

  // check make homogeneous mode space
  Eigen::MatrixXd homogeneous_mode_space =
      make_homogeneous_mode_space(dof_space_0);
  // std::cout << "homogeneous_mode_space: \n"
  //           << homogeneous_mode_space << std::endl;
  EXPECT_EQ(homogeneous_mode_space.rows(), 8);
  EXPECT_EQ(homogeneous_mode_space.cols(), 2);

  // check exclude homogeneous mode space
  DoFSpace dof_space_1 = exclude_homogeneous_mode_space(dof_space_0);
  // std::cout << "excluding homogeneous_mode_space: \n"
  //           << dof_space_1.basis() << std::endl;
  EXPECT_EQ(dof_space_1.basis().rows(), 8);
  EXPECT_EQ(dof_space_1.basis().cols(), 6);

  // check make symmery adapted dof space
  SupercellSymInfo const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  bool calc_wedges = false;
  std::optional<VectorSpaceSymReport> sym_report;
  DoFSpace dof_space_2 = make_symmetry_adapted_dof_space(
      dof_space_1, sym_info, invariant_group, calc_wedges, sym_report);
  std::cout << "excluding homogeneous_mode_space, symmetry adapted: \n"
            << dof_space_2.basis() << std::endl;
  EXPECT_EQ(dof_space_2.basis().rows(), dof_space_1.basis().rows());
  EXPECT_EQ(dof_space_2.basis().cols(), dof_space_1.basis().cols());

  // // check symmetry report
  // jsonParser dof_space_json;
  // to_json(dof_space_2, dof_space_json, "test", config_input, sym_report);
  // std::cout << dof_space_json << std::endl;
}

class RestrictedLocalDoFSpaceTest2 : public testing::Test {
 protected:
  std::shared_ptr<CASM::Structure const> shared_prim;
  std::shared_ptr<CASM::Supercell> shared_supercell;  // conventional unit cell

  static xtal::BasicStructure make_prim();

  RestrictedLocalDoFSpaceTest2()
      : shared_prim(std::make_shared<CASM::Structure const>(make_prim())),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, Eigen::Matrix3l::Identity())) {}
};

xtal::BasicStructure RestrictedLocalDoFSpaceTest2::make_prim() {
  // BCC base structure,
  // - with corner atoms {A, B} allowed to displace in xz
  // - with body-centered atoms {A, B} allowed to displace in yz,
  // -> test case where all atoms may displace,
  //    but rigid translations are only allowed in 1d (z)
  using namespace xtal;

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  SiteDoFSet disp_xz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d0x", "d0z"},          // axes names
      Eigen::MatrixXd({        // basis: allow displacements in xz
                       {1., 0., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  SiteDoFSet disp_yz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d1y", "d1z"},          // axes names
      Eigen::MatrixXd({        // basis: allow displacements in yz
                       {0., 1., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  Lattice lat{Eigen::Vector3d{4.0, 0.0, 0.0}, Eigen::Vector3d{0.0, 4.0, 0.0},
              Eigen::Vector3d{0.0, 0.0, 4.0}};

  BasicStructure struc{lat};
  struc.set_basis(
      {Site{Coordinate{0.0, 0.0, 0.0, lat, FRAC}, {A, B}, {disp_xz}},
       Site{Coordinate{0.5, 0.5, 0.5, lat, FRAC}, {A, B}, {disp_yz}}});
  return struc;
}

TEST_F(RestrictedLocalDoFSpaceTest2, FactorGroupSize) {
  auto const &factor_group = shared_prim->factor_group();

  SymInfoOptions opt{CART};
  brief_description(log(), factor_group, shared_prim->lattice(), opt);

  EXPECT_EQ(factor_group.size(), 16);
}

TEST_F(RestrictedLocalDoFSpaceTest2, PrimSiteDoFSymReps) {
  auto const &factor_group = shared_prim->factor_group();

  // check dimensions of prim disp symrep matrices, should be 2x2
  for (auto const &site_symrep_IDs : shared_prim->site_dof_symrep_IDs()) {
    SymGroupRep const &rep =
        factor_group.representation(site_symrep_IDs.find("disp")->second);
    for (SymOp const &op : factor_group) {
      Eigen::MatrixXd symrep_matrix = *(rep.MatrixXd(op));
      EXPECT_EQ(symrep_matrix.rows(), 2);
      EXPECT_EQ(symrep_matrix.cols(), 2);
    }
  }
}

TEST_F(RestrictedLocalDoFSpaceTest2, CollectiveDoFSymReps) {
  // check dimensions of collective dof symrep matrices:
  // for all sites, should be 4x4

  // Construct the disp DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);
  SupercellSymInfo const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  EXPECT_EQ(invariant_group.size(), 16);

  MasterSymGroup symrep_master_group;
  SymGroupRepID id;
  SymGroupRep const &rep = make_dof_space_symrep(
      dof_space, sym_info, invariant_group, symrep_master_group, id);
  EXPECT_EQ(rep.size(), 16);

  std::cout << "collective_dof_symrep" << std::endl;
  Index i = 0;
  for (SymOp const &op : symrep_master_group) {
    Eigen::MatrixXd matrix_rep = *(rep.MatrixXd(op));
    if (i == 0) {
      EXPECT_TRUE(matrix_rep.isIdentity(TOL));
    }
    std::cout << "---" << std::endl;
    std::cout << "i: " << i << "  op.index(): " << op.index() << std::endl;
    std::cout << matrix_rep << std::endl;
    EXPECT_EQ(matrix_rep.rows(), 4);
    EXPECT_EQ(matrix_rep.cols(), 4);
    i++;
  }
}

TEST_F(RestrictedLocalDoFSpaceTest2, ExcludeHomogeneousModeSpace) {
  // print_local_dof_symreps(shared_supercell->sym_info());

  // Construct the restricted disp DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space_0 = make_dof_space(dof_key, config_input);
  std::cout << "including homogeneous_mode_space: \n"
            << dof_space_0.basis() << std::endl;
  EXPECT_EQ(dof_space_0.basis().rows(), 4);
  EXPECT_EQ(dof_space_0.basis().cols(), 4);

  // check make homogeneous mode space
  Eigen::MatrixXd homogeneous_mode_space =
      make_homogeneous_mode_space(dof_space_0);
  std::cout << "homogeneous_mode_space: \n"
            << homogeneous_mode_space << std::endl;
  EXPECT_EQ(homogeneous_mode_space.rows(), 4);
  EXPECT_EQ(homogeneous_mode_space.cols(), 1);

  // check exclude homogeneous mode space
  DoFSpace dof_space_1 = exclude_homogeneous_mode_space(dof_space_0);
  std::cout << "excluding homogeneous_mode_space: \n"
            << dof_space_1.basis() << std::endl;
  EXPECT_EQ(dof_space_1.basis().rows(), 4);
  EXPECT_EQ(dof_space_1.basis().cols(), 3);

  // check make symmery adapted dof space
  SupercellSymInfo const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  bool calc_wedges = false;
  std::optional<VectorSpaceSymReport> sym_report;
  DoFSpace dof_space_2 = make_symmetry_adapted_dof_space(
      dof_space_1, sym_info, invariant_group, calc_wedges, sym_report);
  std::cout << "excluding homogeneous_mode_space, symmetry adapted: \n"
            << dof_space_2.basis() << std::endl;
  EXPECT_EQ(dof_space_2.basis().rows(), dof_space_1.basis().rows());
  EXPECT_EQ(dof_space_2.basis().cols(), dof_space_1.basis().cols());

  // check symmetry report
  jsonParser dof_space_json;
  to_json(dof_space_2, dof_space_json, "test", config_input, sym_report);
  std::cout << dof_space_json << std::endl;
}

class RestrictedLocalDoFSpaceTest3 : public testing::Test {
 protected:
  std::shared_ptr<CASM::Structure const> shared_prim;
  std::shared_ptr<CASM::Supercell> shared_supercell;  // conventional unit cell

  static xtal::BasicStructure make_prim();

  RestrictedLocalDoFSpaceTest3()
      : shared_prim(std::make_shared<CASM::Structure const>(make_prim())),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, Eigen::Matrix3l::Identity())) {}
};

xtal::BasicStructure RestrictedLocalDoFSpaceTest3::make_prim() {
  // BCC base structure,
  // - with corner atoms {A, B} allowed to displace in xz
  // - with body-centered atoms {C, D} allowed to displace in yz,
  // -> test case where all atoms may displace,
  //    but rigid translations are only allowed in 1d (z),
  //    and the asymmetric unit > 1
  using namespace xtal;

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");
  Molecule C = Molecule::make_atom("C");
  Molecule D = Molecule::make_atom("D");

  SiteDoFSet disp_xz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d0x", "d0z"},          // axes names
      Eigen::MatrixXd({        // basis: allow displacements in xz
                       {1., 0., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  SiteDoFSet disp_yz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d1y", "d1z"},          // axes names
      Eigen::MatrixXd({        // basis: allow displacements in yz
                       {0., 1., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  Lattice lat{Eigen::Vector3d{4.0, 0.0, 0.0}, Eigen::Vector3d{0.0, 4.0, 0.0},
              Eigen::Vector3d{0.0, 0.0, 4.0}};

  BasicStructure struc{lat};
  struc.set_basis(
      {Site{Coordinate{0.0, 0.0, 0.0, lat, FRAC}, {A, B}, {disp_xz}},
       Site{Coordinate{0.5, 0.5, 0.5, lat, FRAC}, {C, D}, {disp_yz}}});
  return struc;
}

TEST_F(RestrictedLocalDoFSpaceTest3, FactorGroupSize) {
  auto const &factor_group = shared_prim->factor_group();

  SymInfoOptions opt{CART};
  brief_description(log(), factor_group, shared_prim->lattice(), opt);

  EXPECT_EQ(factor_group.size(), 8);
}

TEST_F(RestrictedLocalDoFSpaceTest3, PrimSiteDoFSymReps) {
  auto const &factor_group = shared_prim->factor_group();

  // check dimensions of prim disp symrep matrices, should be 2x2
  for (auto const &site_symrep_IDs : shared_prim->site_dof_symrep_IDs()) {
    SymGroupRep const &rep =
        factor_group.representation(site_symrep_IDs.find("disp")->second);
    for (SymOp const &op : factor_group) {
      Eigen::MatrixXd symrep_matrix = *(rep.MatrixXd(op));
      EXPECT_EQ(symrep_matrix.rows(), 2);
      EXPECT_EQ(symrep_matrix.cols(), 2);
    }
  }
}

TEST_F(RestrictedLocalDoFSpaceTest3, CollectiveDoFSymReps) {
  // check dimensions of collective dof symrep matrices:
  // for all sites, should be 4x4

  // Construct the disp DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);
  SupercellSymInfo const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  EXPECT_EQ(invariant_group.size(), 8);

  MasterSymGroup symrep_master_group;
  SymGroupRepID id;
  SymGroupRep const &rep = make_dof_space_symrep(
      dof_space, sym_info, invariant_group, symrep_master_group, id);
  EXPECT_EQ(rep.size(), 8);

  std::cout << "collective_dof_symrep" << std::endl;
  Index i = 0;
  for (SymOp const &op : symrep_master_group) {
    Eigen::MatrixXd matrix_rep = *(rep.MatrixXd(op));
    if (i == 0) {
      EXPECT_TRUE(matrix_rep.isIdentity(TOL));
    }
    std::cout << "---" << std::endl;
    std::cout << "i: " << i << "  op.index(): " << op.index() << std::endl;
    std::cout << matrix_rep << std::endl;
    EXPECT_EQ(matrix_rep.rows(), 4);
    EXPECT_EQ(matrix_rep.cols(), 4);
    i++;
  }
}

TEST_F(RestrictedLocalDoFSpaceTest3, ExcludeHomogeneousModeSpace) {
  // print_local_dof_symreps(shared_supercell->sym_info());

  // Construct the restricted disp DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space_0 = make_dof_space(dof_key, config_input);
  std::cout << "including homogeneous_mode_space: \n"
            << dof_space_0.basis() << std::endl;
  EXPECT_EQ(dof_space_0.basis().rows(), 4);
  EXPECT_EQ(dof_space_0.basis().cols(), 4);

  // check make homogeneous mode space
  Eigen::MatrixXd homogeneous_mode_space =
      make_homogeneous_mode_space(dof_space_0);
  std::cout << "homogeneous_mode_space: \n"
            << homogeneous_mode_space << std::endl;
  EXPECT_EQ(homogeneous_mode_space.rows(), 4);
  EXPECT_EQ(homogeneous_mode_space.cols(), 1);

  // check exclude homogeneous mode space
  DoFSpace dof_space_1 = exclude_homogeneous_mode_space(dof_space_0);
  std::cout << "excluding homogeneous_mode_space: \n"
            << dof_space_1.basis() << std::endl;
  EXPECT_EQ(dof_space_1.basis().rows(), 4);
  EXPECT_EQ(dof_space_1.basis().cols(), 3);

  // check make symmery adapted dof space
  SupercellSymInfo const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  bool calc_wedges = false;
  std::optional<VectorSpaceSymReport> sym_report;
  DoFSpace dof_space_2 = make_symmetry_adapted_dof_space(
      dof_space_1, sym_info, invariant_group, calc_wedges, sym_report);
  std::cout << "excluding homogeneous_mode_space, symmetry adapted: \n"
            << dof_space_2.basis() << std::endl;
  EXPECT_EQ(dof_space_2.basis().rows(), dof_space_1.basis().rows());
  EXPECT_EQ(dof_space_2.basis().cols(), dof_space_1.basis().cols());

  // check symmetry report
  jsonParser dof_space_json;
  to_json(dof_space_2, dof_space_json, "test", config_input, sym_report);
  std::cout << dof_space_json << std::endl;
}

class VariableLocalDoFSpaceTest1 : public testing::Test {
 protected:
  std::shared_ptr<CASM::Structure const> shared_prim;
  std::shared_ptr<CASM::Supercell> shared_supercell;  // conventional unit cell

  static xtal::BasicStructure make_prim();

  VariableLocalDoFSpaceTest1()
      : shared_prim(std::make_shared<CASM::Structure const>(make_prim())),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, Eigen::Matrix3l::Identity())) {}
};

xtal::BasicStructure VariableLocalDoFSpaceTest1::make_prim() {
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
  return struc;
}

TEST_F(VariableLocalDoFSpaceTest1, FactorGroupSize) {
  auto const &factor_group = shared_prim->factor_group();

  // COORD_TYPE mode = FRAC;
  // bool include_va = false;
  // jsonParser prim_json;
  // write_prim(*shared_prim, prim_json, mode, include_va);
  // log() << prim_json << std::endl;
  //
  // SymInfoOptions opt{CART};
  // brief_description(log(), factor_group, shared_prim->lattice(), opt);

  EXPECT_EQ(factor_group.size(), 48);
}

TEST_F(VariableLocalDoFSpaceTest1, PrimSiteDoFSymReps) {
  auto const &factor_group = shared_prim->factor_group();

  // check dimensions of prim disp symrep matrices, should depend on
  // sublattice: 3x3, 2x2, 2x2, 2x2
  Index b = 0;
  for (auto const &site_symrep_IDs : shared_prim->site_dof_symrep_IDs()) {
    SymGroupRep const &rep =
        factor_group.representation(site_symrep_IDs.find("disp")->second);
    for (SymOp const &op : factor_group) {
      Eigen::MatrixXd symrep_matrix = *(rep.MatrixXd(op));

      if (b == 0) {
        EXPECT_EQ(symrep_matrix.rows(), 3);
        EXPECT_EQ(symrep_matrix.cols(), 3);
      } else {
        EXPECT_EQ(symrep_matrix.rows(), 2);
        EXPECT_EQ(symrep_matrix.cols(), 2);
      }
    }
    b++;
  }
}

TEST_F(VariableLocalDoFSpaceTest1, CollectiveDoFSymReps) {
  // check dimensions of collective dof symrep matrices:
  // for all sites, should be 9x9 (9=3+2+2+2)

  // Construct the disp DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);
  SupercellSymInfo const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  EXPECT_EQ(invariant_group.size(), 48);

  MasterSymGroup symrep_master_group;
  SymGroupRepID id;
  SymGroupRep const &rep = make_dof_space_symrep(
      dof_space, sym_info, invariant_group, symrep_master_group, id);
  EXPECT_EQ(rep.size(), 48);

  std::cout << "collective_dof_symrep " << std::endl;
  Index i = 0;
  for (SymOp const &op : symrep_master_group) {
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
    EXPECT_EQ(matrix_rep.rows(), 9);
    EXPECT_EQ(matrix_rep.cols(), 9);
    i++;
  }
}

TEST_F(VariableLocalDoFSpaceTest1, SubsetCollectiveDoFSymReps_1) {
  // check dimensions of collective dof symrep matrices:
  // for sites {1, 2, 3}, should be 6x6

  // Construct the disp DoF space, but only on sites 1, 2, 3
  ConfigEnumInput config_input{*shared_supercell, {1, 2, 3}};
  DoFKey dof_key = "disp";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);
  SupercellSymInfo const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  EXPECT_EQ(invariant_group.size(), 48);

  MasterSymGroup symrep_master_group;
  SymGroupRepID id;
  SymGroupRep const &rep = make_dof_space_symrep(
      dof_space, sym_info, invariant_group, symrep_master_group, id);
  EXPECT_EQ(rep.size(), 48);

  std::cout << "collective_dof_symrep" << std::endl;
  Index i = 0;
  for (SymOp const &op : symrep_master_group) {
    Eigen::MatrixXd matrix_rep = *(rep.MatrixXd(op));
    if (i == 0) {
      EXPECT_TRUE(matrix_rep.isIdentity(TOL));
    }
    // std::cout << "---" << std::endl;
    // std::cout << "i: " << i << "  op.index(): " << op.index() << std::endl;
    // std::cout << "symop: "
    //   << brief_description(op, shared_prim->lattice(), SymInfoOptions{CART})
    //   << std::endl;
    // std::cout << "matrix_rep: \n" << matrix_rep << std::endl;
    EXPECT_EQ(matrix_rep.rows(), 6);
    EXPECT_EQ(matrix_rep.cols(), 6);
    i++;
  }
}

TEST_F(VariableLocalDoFSpaceTest1, SubsetCollectiveDoFSymReps_2) {
  // check dimensions of collective dof symrep matrices:
  // for sites {0, 3}, should be 5x5

  // Construct the disp DoF space, but only on sites 0 and 3
  ConfigEnumInput config_input{*shared_supercell, {0, 3}};
  DoFKey dof_key = "disp";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);
  SupercellSymInfo const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  EXPECT_EQ(invariant_group.size(), 16);

  MasterSymGroup symrep_master_group;
  SymGroupRepID id;
  SymGroupRep const &rep = make_dof_space_symrep(
      dof_space, sym_info, invariant_group, symrep_master_group, id);
  EXPECT_EQ(rep.size(), 16);

  std::cout << "collective_dof_symrep" << std::endl;
  Index i = 0;
  for (SymOp const &op : symrep_master_group) {
    Eigen::MatrixXd matrix_rep = *(rep.MatrixXd(op));
    if (i == 0) {
      EXPECT_TRUE(matrix_rep.isIdentity(TOL));
    }
    // std::cout << "---" << std::endl;
    // std::cout << "i: " << i << "  op.index(): " << op.index() << std::endl;
    // std::cout << "symop: "
    //   << brief_description(op, shared_prim->lattice(), SymInfoOptions{CART})
    //   << std::endl;
    // std::cout << "matrix_rep: \n" << matrix_rep << std::endl;
    EXPECT_EQ(matrix_rep.rows(), 5);
    EXPECT_EQ(matrix_rep.cols(), 5);
    i++;
  }
}

TEST_F(VariableLocalDoFSpaceTest1, ExcludeHomogeneousModeSpace) {
  // In this structure, all sites allow displacements, but no rigid
  // translations are possible

  // print_local_dof_symreps(shared_supercell->sym_info());

  // Construct the restricted disp DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space_0 = make_dof_space(dof_key, config_input);
  std::cout << "including homogeneous_mode_space: \n"
            << dof_space_0.basis() << std::endl;
  EXPECT_EQ(dof_space_0.basis().rows(), 9);
  EXPECT_EQ(dof_space_0.basis().cols(), 9);

  // check make homogeneous mode space
  Eigen::MatrixXd homogeneous_mode_space =
      make_homogeneous_mode_space(dof_space_0);
  std::cout << "homogeneous_mode_space: \n"
            << homogeneous_mode_space << std::endl;
  EXPECT_EQ(homogeneous_mode_space.rows(), 9);
  EXPECT_EQ(homogeneous_mode_space.cols(), 0);

  // check exclude homogeneous mode space
  DoFSpace dof_space_1 = exclude_homogeneous_mode_space(dof_space_0);
  std::cout << "excluding homogeneous_mode_space: \n"
            << dof_space_1.basis() << std::endl;
  EXPECT_EQ(dof_space_1.basis().rows(), 9);
  EXPECT_EQ(dof_space_1.basis().cols(), 9);

  // DoFSpace const &dof_space_1 = dof_space_0;

  // check make symmery adapted dof space
  SupercellSymInfo const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  bool calc_wedges = false;
  std::optional<VectorSpaceSymReport> sym_report;
  DoFSpace dof_space_2 = make_symmetry_adapted_dof_space(
      dof_space_1, sym_info, invariant_group, calc_wedges, sym_report);
  std::cout << "excluding homogeneous_mode_space, symmetry adapted: \n"
            << dof_space_2.basis() << std::endl;
  EXPECT_EQ(dof_space_2.basis().rows(), 9);
  EXPECT_EQ(dof_space_2.basis().cols(), 9);

  // check symmetry report
  jsonParser dof_space_json;
  to_json(dof_space_2, dof_space_json, "test", config_input, sym_report);
  std::cout << dof_space_json << std::endl;
}

class DebugLocalDoFSpaceTest : public testing::Test {
 protected:
  std::shared_ptr<CASM::Structure const> shared_prim;  // must make in test
  std::shared_ptr<CASM::Supercell> shared_supercell;   // must make in test

  static xtal::BasicStructure make_prim();

  DebugLocalDoFSpaceTest() {}

  void check_FactorGroupSize(Index factor_group_size);

  void check_PrimSiteDoFSymReps(
      std::vector<std::pair<Index, Index>> sublat_rep_shape);

  void check_CollectiveDoFSymReps(Index invariant_group_size, Index rep_size,
                                  std::pair<Index, Index> rep_shape);

  void check_SubsetCollectiveDoFSymReps(std::set<Index> sites,
                                        Index invariant_group_size,
                                        Index rep_size,
                                        std::pair<Index, Index> rep_shape);

  void check_SymmetryAdaptedDoFSpace(
      std::pair<Index, Index> initial_dof_space_shape,
      std::pair<Index, Index> symmetry_adapted_dof_space_shape);

  void check_ExcludeHomogeneousModeSpace(
      std::pair<Index, Index> initial_dof_space_shape,
      std::pair<Index, Index> homogeneous_mode_space_shape,
      std::pair<Index, Index> dof_space_shape_excluding_homogeneous_modes,
      std::pair<Index, Index> symmetry_adapted_dof_space_shape);
};

void DebugLocalDoFSpaceTest::check_FactorGroupSize(Index factor_group_size) {
  auto const &factor_group = shared_prim->factor_group();

  COORD_TYPE mode = FRAC;
  bool include_va = false;
  jsonParser prim_json;
  write_prim(*shared_prim, prim_json, mode, include_va);
  log() << prim_json << std::endl;

  std::cout << "prim factor group: " << std::endl;
  SymInfoOptions opt{CART};
  brief_description(log(), factor_group, shared_prim->lattice(), opt);

  EXPECT_EQ(factor_group.size(), factor_group_size);
}

void DebugLocalDoFSpaceTest::check_PrimSiteDoFSymReps(
    std::vector<std::pair<Index, Index>> sublat_rep_shape) {
  auto const &factor_group = shared_prim->factor_group();

  // check dimensions of prim disp symrep matrices, should depend on
  // sublattice: 3x3, 2x2, 2x2, 2x2
  Index b = 0;
  for (auto const &site_symrep_IDs : shared_prim->site_dof_symrep_IDs()) {
    SymGroupRep const &rep =
        factor_group.representation(site_symrep_IDs.find("disp")->second);
    for (SymOp const &op : factor_group) {
      Eigen::MatrixXd symrep_matrix = *(rep.MatrixXd(op));

      EXPECT_EQ(symrep_matrix.rows(), sublat_rep_shape[b].first);
      EXPECT_EQ(symrep_matrix.cols(), sublat_rep_shape[b].second);
    }
    b++;
  }
}

void DebugLocalDoFSpaceTest::check_CollectiveDoFSymReps(
    Index invariant_group_size, Index rep_size,
    std::pair<Index, Index> rep_shape) {
  // check dimensions of collective dof symrep matrices:
  // for all sites, should be 9x9 (9=3+2+2+2)

  std::cout << "supercell factor group: " << std::endl;
  SymInfoOptions opt{CART};
  brief_description(log(), shared_supercell->sym_info().factor_group(),
                    shared_supercell->sym_info().supercell_lattice(), opt);

  // Construct the disp DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);
  SupercellSymInfo const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  EXPECT_EQ(invariant_group.size(), invariant_group_size);

  MasterSymGroup symrep_master_group;
  SymGroupRepID id;
  SymGroupRep const &rep = make_dof_space_symrep(
      dof_space, sym_info, invariant_group, symrep_master_group, id);
  EXPECT_EQ(rep.size(), rep_size);

  std::cout << "collective_dof_symrep " << std::endl;
  Index i = 0;
  for (SymOp const &op : symrep_master_group) {
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
    EXPECT_EQ(matrix_rep.rows(), rep_shape.first);
    EXPECT_EQ(matrix_rep.cols(), rep_shape.second);
    i++;
  }
}

void DebugLocalDoFSpaceTest::check_SubsetCollectiveDoFSymReps(
    std::set<Index> sites, Index invariant_group_size, Index rep_size,
    std::pair<Index, Index> rep_shape) {
  // check dimensions of collective dof symrep matrices:
  // for sites {0, 3}, should be 5x5

  // Construct the disp DoF space, but only on sites 0 and 3
  ConfigEnumInput config_input{*shared_supercell, sites};
  DoFKey dof_key = "disp";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);
  SupercellSymInfo const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  EXPECT_EQ(invariant_group.size(), invariant_group_size);

  MasterSymGroup symrep_master_group;
  SymGroupRepID id;
  SymGroupRep const &rep = make_dof_space_symrep(
      dof_space, sym_info, invariant_group, symrep_master_group, id);
  EXPECT_EQ(rep.size(), rep_size);

  std::cout << "collective_dof_symrep" << std::endl;
  Index i = 0;
  for (SymOp const &op : symrep_master_group) {
    Eigen::MatrixXd matrix_rep = *(rep.MatrixXd(op));
    if (i == 0) {
      EXPECT_TRUE(matrix_rep.isIdentity(TOL));
    }
    // std::cout << "---" << std::endl;
    // std::cout << "i: " << i << "  op.index(): " << op.index() << std::endl;
    // std::cout << "symop: "
    //   << brief_description(op, shared_prim->lattice(), SymInfoOptions{CART})
    //   << std::endl;
    // std::cout << "matrix_rep: \n" << matrix_rep << std::endl;
    EXPECT_EQ(matrix_rep.rows(), rep_shape.first);
    EXPECT_EQ(matrix_rep.cols(), rep_shape.second);
    i++;
  }
}

void DebugLocalDoFSpaceTest::check_SymmetryAdaptedDoFSpace(
    std::pair<Index, Index> initial_dof_space_shape,
    std::pair<Index, Index> symmetry_adapted_dof_space_shape) {
  // In this structure, all sites allow displacements, but no rigid
  // translations are possible

  // print_local_dof_symreps(shared_supercell->sym_info());

  // Construct the restricted disp DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);
  std::cout << "dof space basis: \n" << dof_space.basis() << std::endl;
  EXPECT_EQ(dof_space.basis().rows(), initial_dof_space_shape.first);
  EXPECT_EQ(dof_space.basis().cols(), initial_dof_space_shape.second);

  // check make symmery adapted dof space
  SupercellSymInfo const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  bool calc_wedges = false;
  std::optional<VectorSpaceSymReport> sym_report;
  DoFSpace dof_space_1 = make_symmetry_adapted_dof_space(
      dof_space, sym_info, invariant_group, calc_wedges, sym_report);
  std::cout << "symmetry adapted dof space basis: \n"
            << dof_space_1.basis() << std::endl;
  EXPECT_EQ(dof_space_1.basis().rows(), symmetry_adapted_dof_space_shape.first);
  EXPECT_EQ(dof_space_1.basis().cols(),
            symmetry_adapted_dof_space_shape.second);

  // check symmetry report
  jsonParser dof_space_json;
  to_json(dof_space_1, dof_space_json, "test", config_input, sym_report);
  std::cout << dof_space_json << std::endl;
}

void DebugLocalDoFSpaceTest::check_ExcludeHomogeneousModeSpace(
    std::pair<Index, Index> initial_dof_space_shape,
    std::pair<Index, Index> homogeneous_mode_space_shape,
    std::pair<Index, Index> dof_space_shape_excluding_homogeneous_modes,
    std::pair<Index, Index> symmetry_adapted_dof_space_shape) {
  // In this structure, all sites allow displacements, but no rigid
  // translations are possible

  // print_local_dof_symreps(shared_supercell->sym_info());

  // Construct the restricted disp DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space_0 = make_dof_space(dof_key, config_input);
  std::cout << "including homogeneous_mode_space: \n"
            << dof_space_0.basis() << std::endl;
  EXPECT_EQ(dof_space_0.basis().rows(), initial_dof_space_shape.first);
  EXPECT_EQ(dof_space_0.basis().cols(), initial_dof_space_shape.second);

  // check make homogeneous mode space
  Eigen::MatrixXd homogeneous_mode_space =
      make_homogeneous_mode_space(dof_space_0);
  std::cout << "homogeneous_mode_space: \n"
            << homogeneous_mode_space << std::endl;
  EXPECT_EQ(homogeneous_mode_space.rows(), homogeneous_mode_space_shape.first);
  EXPECT_EQ(homogeneous_mode_space.cols(), homogeneous_mode_space_shape.second);

  // check exclude homogeneous mode space
  DoFSpace dof_space_1 = exclude_homogeneous_mode_space(dof_space_0);
  std::cout << "excluding homogeneous_mode_space: \n"
            << dof_space_1.basis() << std::endl;
  EXPECT_EQ(dof_space_1.basis().rows(),
            dof_space_shape_excluding_homogeneous_modes.first);
  EXPECT_EQ(dof_space_1.basis().cols(),
            dof_space_shape_excluding_homogeneous_modes.second);

  // DoFSpace const &dof_space_1 = dof_space_0;

  // check make symmery adapted dof space
  SupercellSymInfo const &sym_info = shared_supercell->sym_info();
  std::vector<PermuteIterator> invariant_group =
      make_invariant_subgroup(config_input);
  bool calc_wedges = false;
  std::optional<VectorSpaceSymReport> sym_report;
  DoFSpace dof_space_2 = make_symmetry_adapted_dof_space(
      dof_space_1, sym_info, invariant_group, calc_wedges, sym_report);
  std::cout << "excluding homogeneous_mode_space, symmetry adapted: \n"
            << dof_space_2.basis() << std::endl;
  EXPECT_EQ(dof_space_2.basis().rows(), symmetry_adapted_dof_space_shape.first);
  EXPECT_EQ(dof_space_2.basis().cols(),
            symmetry_adapted_dof_space_shape.second);

  // check symmetry report
  jsonParser dof_space_json;
  to_json(dof_space_2, dof_space_json, "test", config_input, sym_report);
  std::cout << dof_space_json << std::endl;
}

TEST_F(DebugLocalDoFSpaceTest, Test1) {  // currently fails
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

  check_FactorGroupSize(48  // Index factor_group_size
  );

  check_PrimSiteDoFSymReps(
      // std::vector<std::pair<Index, Index>> sublat_rep_shape
      {{3, 3}, {2, 2}, {2, 2}, {2, 2}});

  check_CollectiveDoFSymReps(48,     // Index invariant_group_size,
                             48,     // Index rep_size,
                             {9, 9}  // std::pair<Index, Index> rep_shape
  );

  check_SubsetCollectiveDoFSymReps({0},    // std::set<Index> sites,
                                   48,     // Index invariant_group_size,
                                   48,     // Index rep_size,
                                   {3, 3}  // std::pair<Index, Index> rep_shape
  );

  check_ExcludeHomogeneousModeSpace(
      // std::pair<Index, Index> initial_dof_space_shape
      {9, 9},
      // std::pair<Index, Index> homogeneous_mode_space_shape
      {9, 0},
      // std::pair<Index, Index> dof_space_shape_excluding_homogeneous_modes
      {9, 9},
      // std::pair<Index, Index> symmetry_adapted_dof_space_shape
      {9, 9});
}

TEST_F(DebugLocalDoFSpaceTest, Test2) {  // currently fails
  // FCC base structure,
  // - no corner atoms
  // - with face atoms {A, B} allowed to displace in the face plane,
  // -> test case where all atoms may displace, but along different dimensions
  //    so no rigid translations are possible
  using namespace xtal;

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  SiteDoFSet disp_xy{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d0x", "d0y"},          // axes names
      Eigen::MatrixXd({        // basis: allow displacements in xy
                       {1., 0., 0.},
                       {0., 1., 0.}})
          .transpose(),
      {}  // excluded_occs
  };

  SiteDoFSet disp_yz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d1y", "d1z"},          // axes names
      Eigen::MatrixXd({        // basis: allow displacements in yz
                       {0., 1., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  SiteDoFSet disp_xz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d2x", "d2z"},          // axes names
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
      {Site{Coordinate{0.5, 0.5, 0.0, lat, FRAC}, {A, B}, {disp_xy}},
       Site{Coordinate{0.0, 0.5, 0.5, lat, FRAC}, {A, B}, {disp_yz}},
       Site{Coordinate{0.5, 0.0, 0.5, lat, FRAC}, {A, B}, {disp_xz}}});

  shared_prim = std::make_shared<CASM::Structure const>(struc);
  shared_supercell = std::make_shared<CASM::Supercell>(
      shared_prim, Eigen::Matrix3l::Identity());

  check_FactorGroupSize(48  // Index factor_group_size
  );

  check_PrimSiteDoFSymReps(
      // std::vector<std::pair<Index, Index>> sublat_rep_shape
      {{2, 2}, {2, 2}, {2, 2}});

  check_CollectiveDoFSymReps(48,     // Index invariant_group_size,
                             48,     // Index rep_size,
                             {6, 6}  // std::pair<Index, Index> rep_shape
  );

  // check_SubsetCollectiveDoFSymReps(
  //   std::set<Index> sites,
  //   Index invariant_group_size,
  //   Index rep_size,
  //   std::pair<Index, Index> rep_shape);

  check_ExcludeHomogeneousModeSpace(
      // std::pair<Index, Index> initial_dof_space_shape
      {6, 6},
      // std::pair<Index, Index> homogeneous_mode_space_shape
      {6, 0},
      // std::pair<Index, Index> dof_space_shape_excluding_homogeneous_modes
      {6, 6},
      // std::pair<Index, Index> symmetry_adapted_dof_space_shape
      {6, 6});
}

TEST_F(DebugLocalDoFSpaceTest, Test3) {  // passes
  // BCC base structure,
  // - corner atoms {A, B} allowed to displace in xyz
  // - body atoms {A, B} allowed to displace in z,
  // -> test case where all atoms may displace, but along different dimensions
  //    so only 1d (z) rigid translations are possible
  using namespace xtal;

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  SiteDoFSet disp_xyz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d0x", "d0y", "d0z"},   // axes names
      Eigen::MatrixXd({        // basis: allow displacements in xyz
                       {1., 0., 0.},
                       {0., 1., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  SiteDoFSet disp_z{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d1z"},                 // axes names
      Eigen::MatrixXd({        // basis: allow displacements in z
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  Lattice lat{Eigen::Vector3d{4.0, 0.0, 0.0}, Eigen::Vector3d{0.0, 4.0, 0.0},
              Eigen::Vector3d{0.0, 0.0, 4.0}};

  BasicStructure struc{lat};
  struc.set_basis(
      {Site{Coordinate{0.0, 0.0, 0.0, lat, FRAC}, {A, B}, {disp_xyz}},
       Site{Coordinate{0.5, 0.5, 0.5, lat, FRAC}, {A, B}, {disp_z}}});

  shared_prim = std::make_shared<CASM::Structure const>(struc);
  shared_supercell = std::make_shared<CASM::Supercell>(
      shared_prim, Eigen::Matrix3l::Identity());

  check_FactorGroupSize(16  // Index factor_group_size
  );

  check_PrimSiteDoFSymReps(
      // std::vector<std::pair<Index, Index>> sublat_rep_shape
      {{3, 3}, {1, 1}});

  check_CollectiveDoFSymReps(16,     // Index invariant_group_size,
                             16,     // Index rep_size,
                             {4, 4}  // std::pair<Index, Index> rep_shape
  );

  // check_SubsetCollectiveDoFSymReps(
  //   std::set<Index> sites,
  //   Index invariant_group_size,
  //   Index rep_size,
  //   std::pair<Index, Index> rep_shape);

  check_ExcludeHomogeneousModeSpace(
      // std::pair<Index, Index> initial_dof_space_shape
      {4, 4},
      // std::pair<Index, Index> homogeneous_mode_space_shape
      {4, 1},
      // std::pair<Index, Index> dof_space_shape_excluding_homogeneous_modes
      {4, 3},
      // std::pair<Index, Index> symmetry_adapted_dof_space_shape
      {4, 3});
}

TEST_F(DebugLocalDoFSpaceTest, Test4) {  // fails for 2x2x2, passes for 2x1x1
  // BCC base structure, non-primitive supercell
  // - corner atoms {A, B} allowed to displace in xyz
  // - body atoms {A, B} allowed to displace in z,
  // -> test case where all atoms may displace, but along different dimensions
  //    so only 1d (z) rigid translations are possible
  using namespace xtal;

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  SiteDoFSet disp_xyz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d0x", "d0y", "d0z"},   // axes names
      Eigen::MatrixXd({        // basis: allow displacements in xyz
                       {1., 0., 0.},
                       {0., 1., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  SiteDoFSet disp_z{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d1z"},                 // axes names
      Eigen::MatrixXd({        // basis: allow displacements in z
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  Lattice lat{Eigen::Vector3d{4.0, 0.0, 0.0}, Eigen::Vector3d{0.0, 4.0, 0.0},
              Eigen::Vector3d{0.0, 0.0, 4.0}};

  BasicStructure struc{lat};
  struc.set_basis(
      {Site{Coordinate{0.0, 0.0, 0.0, lat, FRAC}, {A, B}, {disp_xyz}},
       Site{Coordinate{0.5, 0.5, 0.5, lat, FRAC}, {A, B}, {disp_z}}});

  shared_prim = std::make_shared<CASM::Structure const>(struc);

  Eigen::Matrix3l T;
  T << 2, 0, 0, 0, 2, 0, 0, 0, 2;

  shared_supercell = std::make_shared<CASM::Supercell>(shared_prim, T);

  Index prim_factor_group_size = 16;
  Index invariant_group_size = 8;
  Index vol = T.determinant();
  Index prim_disp_dof_space_dim = 4;
  Index homogeneous_mode_dim = 1;

  check_FactorGroupSize(prim_factor_group_size  // Index factor_group_size
  );

  check_PrimSiteDoFSymReps(
      // std::vector<std::pair<Index, Index>> sublat_rep_shape
      {{3, 3}, {1, 1}});

  check_CollectiveDoFSymReps(
      invariant_group_size * vol,  // Index invariant_group_size,
      invariant_group_size * vol,  // Index rep_size,
      // std::pair<Index, Index> rep_shape
      {prim_disp_dof_space_dim * vol, prim_disp_dof_space_dim * vol});

  check_SymmetryAdaptedDoFSpace(
      // std::pair<Index, Index> initial_dof_space_shape
      {prim_disp_dof_space_dim * vol, prim_disp_dof_space_dim * vol},
      // std::pair<Index, Index> symmetry_adapted_dof_space_shape
      {prim_disp_dof_space_dim * vol, prim_disp_dof_space_dim * vol});

  check_ExcludeHomogeneousModeSpace(
      // std::pair<Index, Index> initial_dof_space_shape
      {prim_disp_dof_space_dim * vol, prim_disp_dof_space_dim * vol},
      // std::pair<Index, Index> homogeneous_mode_space_shape
      {prim_disp_dof_space_dim * vol, homogeneous_mode_dim},
      // std::pair<Index, Index> dof_space_shape_excluding_homogeneous_modes
      {prim_disp_dof_space_dim * vol,
       prim_disp_dof_space_dim * vol - homogeneous_mode_dim},
      // std::pair<Index, Index> symmetry_adapted_dof_space_shape
      {prim_disp_dof_space_dim * vol,
       prim_disp_dof_space_dim * vol - homogeneous_mode_dim});
}

// fails for transformation_matrix_to_super:
//  0  1 -2
//  0  1  2
//  1 -1  0
TEST_F(DebugLocalDoFSpaceTest, Test5) {
  // FCC, 3d displacements, various supercells
  using namespace xtal;

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  SiteDoFSet disp_xyz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"dx", "dy", "dz"},      // axes names
      Eigen::MatrixXd({        // basis: allow displacements in xyz
                       {1., 0., 0.},
                       {0., 1., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  Lattice lat{Eigen::Vector3d{2.0, 2.0, 0.0}, Eigen::Vector3d{0.0, 2.0, 2.0},
              Eigen::Vector3d{2.0, 0.0, 2.0}};

  BasicStructure struc{lat};
  struc.set_basis(
      {Site{Coordinate{0.0, 0.0, 0.0, lat, FRAC}, {A, B}, {disp_xyz}}});

  shared_prim = std::make_shared<CASM::Structure const>(struc);

  Index prim_factor_group_size = 48;
  Index prim_disp_dof_space_dim = 3;
  Index homogeneous_mode_dim = 3;

  // lattice enumeration
  int begin_volume = 2;
  int end_volume = 5;
  std::string dirs = "abc";
  Eigen::Matrix3i generating_matrix = Eigen::Matrix3i::Identity();
  ScelEnumProps enumeration_params{begin_volume, end_volume, dirs,
                                   generating_matrix};
  CASM::ScelEnumByProps enumerator{shared_prim, enumeration_params};

  // for various supercells:
  for (auto const supercell : enumerator) {
    shared_supercell = std::make_shared<CASM::Supercell>(supercell);
    Eigen::Matrix3l T =
        shared_supercell->sym_info().transformation_matrix_to_super();

    std::cout << "--- begin supercell ---" << std::endl;
    std::cout << "transformation_matrix_to_super:" << std::endl;
    std::cout << T << std::endl;

    Index invariant_group_size =
        shared_supercell->sym_info().factor_group().size();
    Index vol = T.determinant();

    check_FactorGroupSize(prim_factor_group_size  // Index factor_group_size
    );

    check_PrimSiteDoFSymReps(
        // std::vector<std::pair<Index, Index>> sublat_rep_shape
        {{3, 3}});

    check_CollectiveDoFSymReps(
        invariant_group_size * vol,  // Index invariant_group_size,
        invariant_group_size * vol,  // Index rep_size,
        // std::pair<Index, Index> rep_shape
        {prim_disp_dof_space_dim * vol, prim_disp_dof_space_dim * vol});

    check_SymmetryAdaptedDoFSpace(
        // std::pair<Index, Index> initial_dof_space_shape
        {prim_disp_dof_space_dim * vol, prim_disp_dof_space_dim * vol},
        // std::pair<Index, Index> symmetry_adapted_dof_space_shape
        {prim_disp_dof_space_dim * vol, prim_disp_dof_space_dim * vol});

    check_ExcludeHomogeneousModeSpace(
        // std::pair<Index, Index> initial_dof_space_shape
        {prim_disp_dof_space_dim * vol, prim_disp_dof_space_dim * vol},
        // std::pair<Index, Index> homogeneous_mode_space_shape
        {prim_disp_dof_space_dim * vol, homogeneous_mode_dim},
        // std::pair<Index, Index> dof_space_shape_excluding_homogeneous_modes
        {prim_disp_dof_space_dim * vol,
         prim_disp_dof_space_dim * vol - homogeneous_mode_dim},
        // std::pair<Index, Index> symmetry_adapted_dof_space_shape
        {prim_disp_dof_space_dim * vol,
         prim_disp_dof_space_dim * vol - homogeneous_mode_dim});
  }
}
