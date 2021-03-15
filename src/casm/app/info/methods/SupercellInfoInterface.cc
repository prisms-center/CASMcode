#include "casm/app/info/methods/SupercellInfoInterface.hh"

#include "casm/app/ProjectSettings.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/casm_io/dataformatter/DataFormatterTools_impl.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/dataformatter/DatumFormatterAdapter.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
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

/// Data structure for formatting properties of prim and supercell_sym_info
struct SupercellInfoData {
  SupercellInfoData(std::shared_ptr<Structure const> _shared_prim,
                    SupercellSymInfo const &_supercell_sym_info)
      : shared_prim(_shared_prim), supercell_sym_info(_supercell_sym_info) {}

  std::shared_ptr<Structure const> shared_prim;

  SupercellSymInfo const &supercell_sym_info;
};

// use for ValueType=bool, int, double, std::string, jsonParser
template <typename ValueType>
using SupercellInfoFormatter =
    GenericDatumFormatter<ValueType, SupercellInfoData>;

SupercellInfoFormatter<std::string> supercell_name() {
  return SupercellInfoFormatter<std::string>(
      "supercell_name",
      "Unique name given to a supercell, based on the hermite normal form of "
      "the transformation_matrix_to_super and, if not canonical, the index of "
      "the prim factor group operation that transforms the canonical supercell "
      "into this supercell.",
      [](SupercellInfoData const &data) -> std::string {
        return make_supercell_name(data.shared_prim->point_group(),
                                   data.supercell_sym_info.prim_lattice(),
                                   data.supercell_sym_info.supercell_lattice());
      });
}

SupercellInfoFormatter<std::string> canonical_supercell_name() {
  return SupercellInfoFormatter<std::string>(
      "canonical_supercell_name", "Name of the canonical equivalent supercell.",
      [](SupercellInfoData const &data) -> std::string {
        return make_canonical_supercell_name(
            data.shared_prim->point_group(),
            data.supercell_sym_info.prim_lattice(),
            data.supercell_sym_info.supercell_lattice());
      });
}

}  // namespace

namespace adapter {

template <typename ToType, typename FromType>
struct Adapter;

template <>
struct Adapter<SupercellSymInfo, SupercellInfoData> {
  SupercellSymInfo const &operator()(SupercellInfoData const &adaptable) const {
    return adaptable.supercell_sym_info;
  }
};

}  // namespace adapter

namespace {

DataFormatterDictionary<SupercellInfoData> make_supercell_info_dict() {
  DataFormatterDictionary<SupercellInfoData> supercell_info_dict;

  // properties that require prim and supercell_sym_info
  supercell_info_dict.insert(supercell_name(), canonical_supercell_name());

  // properties that only require supercell_sym_info
  auto sym_info_dict = make_dictionary<SupercellSymInfo>();
  for (auto it = sym_info_dict.begin(); it != sym_info_dict.end(); ++it) {
    if (it->type() == DatumFormatterClass::Property) {
      supercell_info_dict.insert(
          make_datum_formatter_adapter<SupercellInfoData, SupercellSymInfo>(
              *it));
    }
  }

  return supercell_info_dict;
}

}  // namespace

std::string SupercellInfoInterface::desc() const {
  std::string description =
      "Get information about a supercell. The supercell is specified by   \n"
      "the prim and one of the following (else the primitive cell is      \n"
      "used):                                                             \n"
      "- transformation_matrix_to_super                                   \n"
      "- supercell_lattice_vectors                                        \n"
      "- supercell_lattice_column_matrix                                  \n"
      "- supercell_name                                                   \n\n";

  std::string custom_options =
      "  prim: JSON object (optional, default=prim of current project)    \n"
      "    See `casm format --prim` for details on the prim format.       \n\n"

      "  transformation_matrix_to_super: 3x3 array of integer (optional)  \n"
      "    Transformation matrix T, defining the supercell lattice vectors\n"
      "    S, in terms of the prim lattice vectors, P: `S = P * T`, where \n"
      "    S and P are column vector matrices.                            \n\n"

      "  supercell_lattice_vectors: 3x3 array of integer (optional)       \n"
      "    Supercell lattice vectors, as a row vector matrix.             \n\n"

      "  supercell_lattice_column_matrix: 3x3 array of integer (optional) \n"
      "    Supercell lattice vectors, as a column vector matrix.          \n\n"

      "  supercell_name: string (optional)                                \n"
      "    Unique name given to a supercell, based on the hermite normal  \n"
      "    form, of the transformation_matrix_to_super and, if not        \n"
      "    canonical, the index of the prim factor group operation that   \n"
      "    transforms the canonical supercell into this supercell.        \n\n"

      "  properties: array of string (optional, default=[])               \n"
      "    An array of strings specifying which supercell properties to   \n"
      "    output. The default value of an empty array will print all     \n"
      "    properties. The allowed options are:                           \n\n";

  std::stringstream ss;
  auto dict = make_supercell_info_dict();
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
  if (parser.self.contains("transformation_matrix_to_super")) {
    parser.optional(T, "transformation_matrix_to_super");

    // or read "supercell_lattice_vectors"
  } else if (parser.self.contains("supercell_lattice_vectors")) {
    Eigen::Matrix3d L_transpose;
    parser.optional(L_transpose, "supercell_lattice_vectors");
    Lattice super_lattice{L_transpose.transpose()};
    T = make_transformation_matrix_to_super(shared_prim->lattice(),
                                            super_lattice, TOL);

    // or read "supercell_lattice_column_matrix"
  } else if (parser.self.contains("supercell_lattice_column_matrix")) {
    Eigen::Matrix3d L;
    parser.optional(L, "supercell_lattice_column_matrix");
    Lattice super_lattice{L};
    T = make_transformation_matrix_to_super(shared_prim->lattice(),
                                            super_lattice, TOL);

    // or read "supercell_name"
  } else if (parser.self.contains("supercell_name")) {
    std::string supercell_name;
    parser.optional(supercell_name, "supercell_name");
    xtal::Superlattice superlattice = make_superlattice_from_supercell_name(
        shared_prim->factor_group(), shared_prim->lattice(), supercell_name);
    T = superlattice.transformation_matrix_to_super();

    // else use Identity (prim cell)
  } else {
    T = Eigen::Matrix3l::Identity();
  }

  // read "properties"
  std::vector<std::string> properties;
  parser.optional(properties, "properties");
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  // construct SupercellSymInfo
  Lattice super_lattice = make_superlattice(shared_prim->lattice(), T);
  SupercellSymInfo supercell_sym_info =
      make_supercell_sym_info(*shared_prim, super_lattice);

  DataFormatterDictionary<SupercellInfoData> supercell_info_dict =
      make_supercell_info_dict();

  // output all properties if empty
  if (properties.empty()) {
    auto it = supercell_info_dict.begin();
    for (; it != supercell_info_dict.end(); ++it) {
      if (it->type() == DatumFormatterClass::Property) {
        properties.push_back(it->name());
      }
    }
  }

  // format
  auto formatter = supercell_info_dict.parse(properties);
  jsonParser json;
  SupercellInfoData data{shared_prim, supercell_sym_info};
  formatter.to_json(data, json);
  log << json << std::endl;
}

}  // namespace CASM
