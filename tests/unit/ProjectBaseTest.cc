#include "ProjectBaseTest.hh"

#include "Common.hh"
#include "autotools.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectBuilder.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/basis_set/BasisFunctionSpecs.hh"
#include "casm/basis_set/OccupationDoFTraits.hh"
#include "casm/clex/ClexBasisSpecs.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex_impl.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/SupercellSymInfo.hh"

namespace test {

ProjectBaseTest::ProjectBaseTest(xtal::BasicStructure const &basic_structure,
                                 std::string title,
                                 jsonParser const &_basis_set_specs_json,
                                 fs::path projpath)
    : tmp_dir(),
      shared_prim(std::make_shared<Structure const>(basic_structure)),
      project_settings_ptr(notstd::make_unique<ProjectSettings>(
          make_default_project_settings(*shared_prim, title, projpath))),
      basis_set_name("default"),
      basis_set_specs_json(_basis_set_specs_json) {
  // for tests, need to set compiler search paths to the CASM repository
  project_settings_ptr->set_casm_libdir(autotools::abs_libdir());
  project_settings_ptr->set_casm_includedir(autotools::abs_includedir());
  if (!fs::exists(projpath)) {
    build_project(*project_settings_ptr, *shared_prim);
  }
  this->write_bspecs_json();

  primclex_ptr = std::make_unique<PrimClex>(project_settings_ptr->root_dir());
}

ProjectBaseTest::ProjectBaseTest(xtal::BasicStructure const &basic_structure,
                                 std::string title,
                                 jsonParser const &_basis_set_specs_json)
    : tmp_dir(),
      shared_prim(std::make_shared<Structure const>(basic_structure)),
      project_settings_ptr(notstd::make_unique<ProjectSettings>(
          make_default_project_settings(*shared_prim, title, tmp_dir.path()))),
      basis_set_name("default"),
      basis_set_specs_json(_basis_set_specs_json) {
  // for tests, need to set compiler search paths to the CASM repository
  project_settings_ptr->set_casm_libdir(autotools::abs_libdir());
  project_settings_ptr->set_casm_includedir(autotools::abs_includedir());
  build_project(*project_settings_ptr, *shared_prim);

  this->write_bspecs_json();

  primclex_ptr = std::make_unique<PrimClex>(project_settings_ptr->root_dir());
}

ProjectBaseTest::~ProjectBaseTest() {}

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
      CASM::make_clexulator(settings, basis_set_name, prim_neighbor_list);
}

int &ProjectBaseTest::occ(Configuration &config,
                          xtal::UnitCellCoord unitcellcoord) {
  auto const &scel_sym_info = config.supercell().sym_info();
  Index l = scel_sym_info.unitcellcoord_index_converter()(unitcellcoord);
  return config.configdof().occ(l);
}

std::vector<IntegralCluster> ProjectBaseTest::make_phenom_orbit() {
  auto shared_prim = primclex_ptr->shared_prim();
  ClexBasisSpecs const &basis_set_specs =
      primclex_ptr->basis_set_specs(basis_set_name);
  ClusterSpecs const &cluster_specs = *basis_set_specs.cluster_specs;
  if (cluster_specs.periodicity_type() != CLUSTER_PERIODICITY_TYPE::LOCAL) {
    throw std::runtime_error(
        "Error in ProjectBaseTest::make_phenom_orbit: is only valid for local "
        "cluster expansions.");
  }
  IntegralCluster prototype_phenom = cluster_specs.get_phenomenal_cluster();
  std::vector<SymOp> equivalents_generating_ops =
      make_equivalents_generating_ops(primclex_ptr->shared_prim(),
                                      prototype_phenom,
                                      cluster_specs.get_generating_group());
  std::vector<IntegralCluster> phenom_orbit;
  for (auto const &op : equivalents_generating_ops) {
    phenom_orbit.push_back(sym::copy_apply(op, prototype_phenom, *shared_prim));
  }
  return phenom_orbit;
}

}  // namespace test
