#include "casm/app/info/methods/PrimInfoInterface.hh"

#include "casm/app/ProjectSettings.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/casm_io/dataformatter/DataFormatterTools_impl.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/crystallography/io/data/SharedPrim_data_io.hh"

namespace CASM {

namespace {

std::shared_ptr<Structure const> open_shared_prim(fs::path root) {
  ProjectSettings settings = open_project_settings(root);
  return std::make_shared<Structure const>(
      read_prim(settings.dir().prim(), settings.crystallography_tol()));
}

}  // namespace

std::string PrimInfoInterface::desc() const {
  std::string description = "Get information about a prim structure.\n\n";

  std::string custom_options =
      "  prim: JSON object (optional, default=prim of current project) \n"
      "    See `casm format --prim` for details on the prim format.\n\n"

      "  properties: array of string                                      \n"
      "    An array of strings specifying which prim properties to output.\n"
      "    The allowed options are:                                       \n\n";

  std::stringstream ss;
  auto dict = make_dictionary<std::shared_ptr<Structure const>>();
  dict.print_help(ss, DatumFormatterClass::Property);

  return name() + ": \n\n" + description + custom_options + ss.str();
}

std::string PrimInfoInterface::name() const { return "PrimInfo"; }

/// Run `prim` info method
void PrimInfoInterface::run(jsonParser const &json_options,
                            PrimClex const *primclex) const {
  Log &log = CASM::log();

  ParentInputParser parser{json_options};
  std::runtime_error error_if_invalid{"Error reading PrimInfo input"};

  std::shared_ptr<Structure const> shared_prim;
  if (parser.self.contains("prim")) {
    // prim provided in input
    xtal::BasicStructure basic_structure;
    parser.optional<xtal::BasicStructure>(basic_structure, "prim", TOL);
    if (parser.valid()) {
      shared_prim = std::make_shared<Structure const>(basic_structure);
    }
  } else if (primclex != nullptr) {
    // if project provided via api
    shared_prim = primclex->shared_prim();
  } else {
    // if project contains current working directory
    fs::path root = find_casmroot(fs::current_path());
    if (!root.empty()) {
      try {
        shared_prim = open_shared_prim(root);
      } catch (std::exception &e) {
        parser.insert_error("prim", e.what());
      }
    } else {
      std::stringstream msg;
      msg << "Error in PrimInfo: No \"prim\" in input and no project provided "
             "or found.";
      parser.insert_error("prim", msg.str());
    }
  }

  std::vector<std::string> properties;
  parser.require(properties, "properties");
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  auto dict = make_dictionary<std::shared_ptr<Structure const>>();
  auto formatter = dict.parse(properties);

  jsonParser json;
  formatter.to_json(shared_prim, json);
  log << json << std::endl;
}

}  // namespace CASM
