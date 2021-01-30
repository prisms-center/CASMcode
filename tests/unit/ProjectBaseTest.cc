#include "ProjectBaseTest.hh"

#include "Common.hh"
#include "autotools.hh"
#include "casm/app/ProjectBuilder.hh"
#include "casm/basis_set/BasisFunctionSpecs.hh"
#include "casm/basis_set/OccupationDoFTraits.hh"
#include "casm/clex/ClexBasisSpecs.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "casm/crystallography/Structure.hh"

namespace test {

ProjectBaseTest::ProjectBaseTest(xtal::BasicStructure const &basic_structure,
                                 std::string title,
                                 jsonParser const &_basis_set_specs_json)
    : shared_prim(std::make_shared<Structure const>(basic_structure)),
      project_settings_ptr(
          notstd::make_unique<ProjectSettings>(make_default_project_settings(
              *shared_prim, title,
              test::proj_dir(autotools::abs_srcdir() +
                             "/tests/unit/test_projects/" + title)))),
      basis_set_name("default"),
      basis_set_specs_json(_basis_set_specs_json) {
  // for tests, need to set compiler search paths to the CASM repository
  project_settings_ptr->set_casm_libdir(autotools::abs_libdir());
  project_settings_ptr->set_casm_includedir(autotools::abs_includedir());
  build_project(*project_settings_ptr, *shared_prim);

  this->write_bspecs_json();

  primclex_ptr = std::make_unique<PrimClex>(project_settings_ptr->root_dir());

  this->write_basis_set_data();
  this->make_clexulator();
}

void ProjectBaseTest::write_bspecs_json() {
  // Write bspecs.json file
  DirectoryStructure const &dir = project_settings_ptr->dir();
  fs::path basis_set_specs_path = dir.bspecs(basis_set_name);
  fs::ofstream file{basis_set_specs_path};
  basis_set_specs_json.print(file);
  file.close();
}

void ProjectBaseTest::write_basis_set_data() {
  // Write clust.json, basis.json, clexulator source file
  auto shared_prim = primclex_ptr->shared_prim();
  auto const &settings = primclex_ptr->settings();
  auto const &basis_set_specs = primclex_ptr->basis_set_specs(basis_set_name);
  auto &prim_neighbor_list = primclex_ptr->nlist();
  CASM::write_basis_set_data(shared_prim, settings, basis_set_name,
                             basis_set_specs, prim_neighbor_list);
}

void ProjectBaseTest::make_clexulator() {
  // Compile clexulator
  auto const &settings = primclex_ptr->settings();
  auto &prim_neighbor_list = primclex_ptr->nlist();
  Clexulator clexulator =
      make_clexulator(settings, basis_set_name, prim_neighbor_list);
}

}  // namespace test
