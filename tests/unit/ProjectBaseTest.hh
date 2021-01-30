#ifndef CASM_unit_ProjectBaseTest
#define CASM_unit_ProjectBaseTest

#include <boost/filesystem.hpp>
#include <string>

#include "casm/global/definitions.hh"

namespace CASM {
namespace xtal {
class BasicStructure;
}
class PrimClex;
class ProjectSettings;
class Structure;
class jsonParser;
}  // namespace CASM

namespace test {

// Creates a project with Clexulator ready to go
class ProjectBaseTest : public testing::Test {
 protected:
  ProjectBaseTest(xtal::BasicStructure const &basic_structure,
                  std::string title, jsonParser const &_basis_set_specs_json);

  void write_bspecs_json();

  void write_basis_set_data();

  void make_clexulator();

  std::shared_ptr<Structure const> shared_prim;
  fs::path root_dir;
  std::unique_ptr<ProjectSettings> project_settings_ptr;
  std::string basis_set_name;
  jsonParser basis_set_specs_json;
  std::unique_ptr<PrimClex> primclex_ptr;
};

}  // namespace test

#endif
