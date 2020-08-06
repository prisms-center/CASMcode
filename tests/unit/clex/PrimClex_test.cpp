#include "gtest/gtest.h"
#include "autotools.hh"

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "clex/TestClexBasisSpecs.hh"

#include "casm/app/ProjectBuilder.hh"
#include "casm/basis_set/BasisFunctionSpecs.hh"
#include "casm/basis_set/OccupationDoFTraits.hh"
#include "casm/clex/ClexBasisSpecs.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "casm/crystallography/Structure.hh"

// PrimClex without directory for saving project data
class FCCTernaryPrimClexTest : public testing::Test {
protected:

  void SetUp() override {

    std::string title = "FCCTernaryPrimClexTest";

    shared_prim = std::make_shared<Structure const>(test::FCC_ternary_prim());
    project_settings_ptr = notstd::make_unique<ProjectSettings>(
                             make_default_project_settings(*shared_prim, title));
  }

  std::shared_ptr<Structure const> shared_prim;
  std::unique_ptr<ProjectSettings> project_settings_ptr;
};

// PrimClex with directory for saving project data
class FCCTernaryProjectTest : public testing::Test {
protected:

  void SetUp() override {

    std::string title = "FCCTernaryProjectTest";

    // Creates a new project directory, appending ".(#)" to ensure it is a new project
    fs::path tmp_root_dir = autotools::abs_srcdir() + "/tests/unit/test_projects/" + title;
    fs::path root_dir = test::proj_dir(tmp_root_dir);

    shared_prim = std::make_shared<Structure const>(test::FCC_ternary_prim());
    project_settings_ptr = notstd::make_unique<ProjectSettings>(
                             make_default_project_settings(*shared_prim, title, root_dir));

    // for tests, need to set compiler search paths to the CASM repository
    project_settings_ptr->set_casm_libdir(autotools::abs_libdir());
    project_settings_ptr->set_casm_includedir(autotools::abs_includedir());

    build_project(*project_settings_ptr, *shared_prim);
  }

  std::shared_ptr<Structure const> shared_prim;
  std::unique_ptr<ProjectSettings> project_settings_ptr;
};

TEST_F(FCCTernaryPrimClexTest, Basics) {

  EXPECT_EQ(shared_prim->basis().size(), 1);

  // Construct PrimClex from project settings & prim
  PrimClex primclex {*project_settings_ptr, shared_prim, log()};
  EXPECT_EQ(primclex.prim().basis().size(), 1);

}

TEST_F(FCCTernaryProjectTest, ReadNewProject) {

  // Construct PrimClex from existing project directory
  PrimClex primclex {project_settings_ptr->root_dir(), log()};
  EXPECT_EQ(primclex.prim().basis().size(), 1);

}

TEST_F(FCCTernaryProjectTest, ReadClexBasisSpecs) {

  // Write bspecs.json file
  DirectoryStructure const &dir = project_settings_ptr->dir();
  std::string basis_set_name = "default";
  fs::path basis_set_specs_path = dir.bspecs(basis_set_name);
  fs::ofstream file {basis_set_specs_path};
  jsonParser basis_set_specs_json = jsonParser::parse(test::FCC_ternary_clex_basis_specs_str_ex0());
  basis_set_specs_json.print(file);
  file.close();

  // Construct PrimClex from existing project directory
  PrimClex primclex {project_settings_ptr->root_dir(), log()};
  EXPECT_EQ(primclex.prim().basis().size(), 1);

  // Assert bspecs.json file exists and is read correctly
  EXPECT_TRUE(primclex.has_basis_set_specs(basis_set_name));

  ClexBasisSpecs const &basis_set_specs = primclex.basis_set_specs(basis_set_name);

  BasisFunctionSpecs const &basis_function_specs = basis_set_specs.basis_function_specs;
  EXPECT_EQ(basis_function_specs.dof_keys.size(), 1);
  EXPECT_EQ(basis_function_specs.dof_keys[0], "occ");

  DoF_impl::OccupationDoFSpecs const &occ_specs = get<DoF_impl::OccupationDoFSpecs>("occ", basis_function_specs);
  EXPECT_EQ(occ_specs.site_basis_function_type, DoF_impl::SITE_BASIS_FUNCTION_TYPE::OCCUPATION);
  EXPECT_EQ(occ_specs.sublat_composition.size(), 0);

  ClusterSpecs const &cluster_specs = *basis_set_specs.cluster_specs;
  EXPECT_EQ(cluster_specs.name(), "periodic_max_length");
  EXPECT_EQ(cluster_specs.periodicity_type(), CLUSTER_PERIODICITY_TYPE::PRIM_PERIODIC);

}

TEST_F(FCCTernaryProjectTest, ClexBasisSet) {

  // Write bspecs.json file
  DirectoryStructure const &dir = project_settings_ptr->dir();
  std::string basis_set_name = "default";
  fs::path basis_set_specs_path = dir.bspecs(basis_set_name);
  fs::ofstream file {basis_set_specs_path};
  jsonParser basis_set_specs_json = jsonParser::parse(test::FCC_ternary_clex_basis_specs_str_ex0());
  basis_set_specs_json.print(file);
  file.close();

  // Construct PrimClex from existing project directory
  PrimClex primclex {project_settings_ptr->root_dir(), log()};
  EXPECT_EQ(primclex.prim().basis().size(), 1);

  // Assert bspecs.json file exists and is read correctly
  EXPECT_TRUE(primclex.has_basis_set_specs(basis_set_name));

  // Write clust.json, basis.json, clexulator source file
  auto shared_prim = primclex.shared_prim();
  auto const &settings = primclex.settings();
  auto const &basis_set_specs = primclex.basis_set_specs(basis_set_name);
  auto &prim_neighbor_list = primclex.nlist();
  write_basis_set_data(shared_prim, settings, basis_set_name, basis_set_specs, prim_neighbor_list);

  // Check that files were created
  auto project_name = settings.project_name();
  EXPECT_TRUE(fs::exists(dir.clust(basis_set_name)));
  EXPECT_TRUE(fs::exists(dir.basis(basis_set_name)));
  EXPECT_TRUE(fs::exists(dir.clexulator_src(project_name, basis_set_name)));

  // Compile clexulator
  Clexulator clexulator = make_clexulator(settings, basis_set_name, prim_neighbor_list);
  EXPECT_TRUE(fs::exists(dir.clexulator_o(project_name, basis_set_name)));
  EXPECT_TRUE(fs::exists(dir.clexulator_so(project_name, basis_set_name)));

  EXPECT_EQ(clexulator.name(), "FCCTernaryProjectTest_Clexulator");
  EXPECT_EQ(clexulator.nlist_size(), 176);
  EXPECT_EQ(clexulator.corr_size(), 75);
  EXPECT_EQ(clexulator.neighborhood().size(), 75);
}
