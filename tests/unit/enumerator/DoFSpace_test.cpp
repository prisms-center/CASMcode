#include "gtest/gtest.h"
#include "autotools.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/enumerator/ConfigEnumInput_impl.hh"
#include "casm/enumerator/DoFSpace_impl.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"
#include "casm/symmetry/io/json/SymRepTools.hh"
#include "crystallography/TestStructures.hh"

using namespace CASM;
using namespace test;

namespace {
  Eigen::Matrix3l _fcc_conventional_transf_mat() {
    Eigen::Matrix3l transf_mat;
    transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    return transf_mat;
  }
}

class DoFSpaceTest : public testing::Test {
protected:

  std::shared_ptr<CASM::Structure const> shared_prim;
  std::shared_ptr<CASM::Supercell> shared_supercell;  // conventional unit cell

  DoFSpaceTest():
    shared_prim(std::make_shared<CASM::Structure const>(test::FCC_ternary_strain_disp_prim())),
    shared_supercell(std::make_shared<CASM::Supercell>(shared_prim, _fcc_conventional_transf_mat())) {}

};

TEST_F(DoFSpaceTest, ConstructorTest1) {
  // Construct the GLstrain DoF space.
  ConfigEnumInput config_input {*shared_supercell};
  DoFKey dof_key = "GLstrain";
  DoFSpace dof_space {config_input, dof_key};
  EXPECT_EQ(dof_space.dof_subspace.rows(), 6);
}

TEST_F(DoFSpaceTest, ConstructorTest2) {
  // Construct the disp DoF space.
  ConfigEnumInput config_input {*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space {config_input, dof_key};
  EXPECT_EQ(dof_space.dof_subspace.rows(), 4 * 3);
}

TEST_F(DoFSpaceTest, JsonIOTest1) {

  // Construct the GLstrain DoF space.
  ConfigEnumInput config_input {*shared_supercell};
  DoFKey dof_key = "GLstrain";
  DoFSpace dof_space {config_input, dof_key};

  jsonParser dof_space_json;
  to_json(dof_space, dof_space_json, shared_supercell->name() + "/0");

  // // Uncomment to print dof_space:
  // log() << dof_space_json << std::endl;

  // a couple tests to check JSON was constructed successfully, but does not check each value
  EXPECT_EQ(dof_space_json.size(), 4);
  EXPECT_EQ(dof_space_json.contains("dof"), true);
  EXPECT_EQ(dof_space_json.contains("initial_configuration"), true);
  EXPECT_EQ(dof_space_json.contains("glossary"), true);
  EXPECT_EQ(dof_space_json.contains("dof_subspace"), true);
}

TEST_F(DoFSpaceTest, JsonIOTest2) {

  // Construct the disp DoF space.
  ConfigEnumInput config_input {*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space {config_input, dof_key};

  jsonParser dof_space_json;
  to_json(dof_space, dof_space_json, shared_supercell->name() + "/0");

  // // Uncomment to print dof_space:
  // log() << dof_space_json << std::endl;

  // a couple tests to check JSON was constructed successfully, but does not check each value
  EXPECT_EQ(dof_space_json.size(), 4);
  EXPECT_EQ(dof_space_json.contains("dof"), true);
  EXPECT_EQ(dof_space_json.contains("initial_configuration"), true);
  EXPECT_EQ(dof_space_json.contains("glossary"), true);
  EXPECT_EQ(dof_space_json.contains("dof_subspace"), true);
}

TEST_F(DoFSpaceTest, VectorSpaceSymReportTest1) {

  // Construct the GLstrain DoF space.
  ConfigEnumInput config_input {*shared_supercell};
  DoFKey dof_key = "GLstrain";
  DoFSpace dof_space {config_input, dof_key};

  // Construct the VectorSpaceSymReport (calc_wedges==false)
  std::vector<PermuteIterator> invariant_group = make_invariant_subgroup(config_input);
  bool calc_wedges = false;
  VectorSpaceSymReport report = vector_space_sym_report(dof_space,
                                                        invariant_group.begin(),
                                                        invariant_group.end(),
                                                        calc_wedges);

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
  ConfigEnumInput config_input {*shared_supercell};
  DoFKey dof_key = "GLstrain";
  DoFSpace dof_space {config_input, dof_key};

  // Construct the VectorSpaceSymReport (calc_wedges==true)
  std::vector<PermuteIterator> invariant_group = make_invariant_subgroup(config_input);
  bool calc_wedges = true;
  VectorSpaceSymReport report = vector_space_sym_report(dof_space,
                                                        invariant_group.begin(),
                                                        invariant_group.end(),
                                                        calc_wedges);
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
  ConfigEnumInput config_input {*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space {config_input, dof_key};

  // Construct the VectorSpaceSymReport (calc_wedges==false)
  std::vector<PermuteIterator> invariant_group = make_invariant_subgroup(config_input);
  bool calc_wedges = false;
  VectorSpaceSymReport report = vector_space_sym_report(dof_space,
                                                        invariant_group.begin(),
                                                        invariant_group.end(),
                                                        calc_wedges);

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
  ConfigEnumInput config_input {*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space {config_input, dof_key};

  // Construct the VectorSpaceSymReport (calc_wedges==true)
  std::vector<PermuteIterator> invariant_group = make_invariant_subgroup(config_input);
  bool calc_wedges = true;
  VectorSpaceSymReport report = vector_space_sym_report(dof_space,
                                                        invariant_group.begin(),
                                                        invariant_group.end(),
                                                        calc_wedges);
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
