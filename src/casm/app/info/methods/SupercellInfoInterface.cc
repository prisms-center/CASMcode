#include "casm/app/info/methods/SupercellInfoInterface.hh"

#include "casm/app/ProjectSettings.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/casm_io/dataformatter/DataFormatterTools_impl.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/symmetry/SupercellSymInfo.hh"
#include "casm/symmetry/io/data/SupercellSymInfo_data_io.hh"

namespace CASM {

namespace {

std::shared_ptr<Structure const> open_shared_prim(fs::path root) {
  ProjectSettings settings = open_project_settings(root);
  return std::make_shared<Structure const>(
      read_prim(settings.dir().prim(), settings.crystallography_tol()));
}

}  // namespace

std::string SupercellInfoInterface::desc() const {
  std::string description = "Get information about a supercell.           \n\n";

  std::string custom_options =
      "  prim: JSON object (optional, default=prim of current project)    \n"
      "    See `casm format --prim` for details on the prim format.       \n\n"

      "  transformation_matrix_to_super: 3x3 array of integer             \n"
      "    Transformation matrix T, defining the supercell lattice vectors\n"
      "    S, in terms of the prim lattice vectors, P: `S = T * P`, where \n"
      "    S and P are column vector matrices.                            \n\n"

      "  properties: array of string (optional, default=[])               \n"
      "    An array of strings specifying which prim properties to output.\n"
      "    The default value of an empty array will print all properties. \n"
      "    The allowed options are:                                       \n\n";

  std::stringstream ss;
  auto dict = make_dictionary<SupercellSymInfo>();
  dict.print_help(ss, DatumFormatterClass::Property);

  return name() + ": \n\n" + description + custom_options + ss.str();
}

std::string SupercellInfoInterface::name() const { return "SupercellInfo"; }

/// Run `prim` info method
void SupercellInfoInterface::run(jsonParser const &json_options,
                                 PrimClex const *primclex) const {
  Log &log = CASM::log();

  ParentInputParser parser{json_options};
  std::runtime_error error_if_invalid{"Error reading SupercellInfo input"};

  // read "prim"
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
      msg << "Error in SupercellInfo: No \"prim\" in input and no project "
             "provided "
             "or found.";
      parser.insert_error("prim", msg.str());
    }
  }

  // read "transformation_matrix_to_super"
  Eigen::Matrix3l T;
  parser.require(T, "transformation_matrix_to_super");

  // read "properties"
  std::vector<std::string> properties;
  parser.optional(properties, "properties");
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  // construct SupercellSymInfo
  Lattice super_lattice = make_superlattice(shared_prim->lattice(), T);
  SupercellSymInfo supercell_sym_info =
      make_supercell_sym_info(*shared_prim, super_lattice);

  auto dict = make_dictionary<SupercellSymInfo>();
  if (properties.empty()) {
    for (auto it = dict.begin(); it != dict.end(); ++it) {
      if (it->type() == DatumFormatterClass::Property) {
        properties.push_back(it->name());
      }
    }
  }
  auto formatter = dict.parse(properties);

  jsonParser json;
  formatter.to_json(supercell_sym_info, json);
  log << json;
}

}  // namespace CASM
