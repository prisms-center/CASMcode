#include "casm/app/info/methods/ConfigurationInfoInterface.hh"

#include "casm/app/ProjectSettings.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/casm_io/dataformatter/DataFormatterTools_impl.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/dataformatter/DatumFormatterAdapter.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/FillSupercell.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/io/json/Configuration_json_io.hh"
#include "casm/clex/io/stream/Configuration_stream_io.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"

namespace CASM {

namespace {

std::shared_ptr<Structure const> open_shared_prim(fs::path root) {
  ProjectSettings settings = open_project_settings(root);
  return std::make_shared<Structure const>(
      read_prim(settings.dir().prim(), settings.crystallography_tol()));
}

/// Data structure for formatting properties of prim and supercell_sym_info
struct ConfigurationInfoData {
  ConfigurationInfoData(PrimClex const *_primclex, fs::path _root,
                        std::shared_ptr<Structure const> _shared_prim,
                        Configuration const &_configuration,
                        jsonParser const &_params_json)
      : shared_prim(_shared_prim),
        configuration(_configuration),
        params_json(_params_json),
        m_primclex(_primclex),
        m_root(_root) {}

  std::shared_ptr<Structure const> shared_prim;

  Configuration const &configuration;

  jsonParser const &params_json;

  PrimClex const &primclex() const {
    if (m_primclex) {
      return *m_primclex;
    } else if (m_root.empty() && !m_uniq_primclex) {
      throw std::runtime_error(
          "Requested project information, but no project exists");
    } else if (!m_uniq_primclex) {
      m_uniq_primclex.reset(new PrimClex(m_root));
    }
    return *m_uniq_primclex;
  }

 private:
  /// May be nullptr
  PrimClex const *m_primclex;
  /// May be empty
  fs::path m_root;
  mutable std::unique_ptr<PrimClex> m_uniq_primclex;
};

// use for ValueType=bool, int, double, std::string, jsonParser
template <typename ValueType>
using ConfigurationInfoFormatter =
    GenericDatumFormatter<ValueType, ConfigurationInfoData>;

typedef Generic1DDatumFormatter<Eigen::VectorXd, ConfigurationInfoData>
    VectorXdConfigurationInfoFormatter;

typedef Generic1DDatumFormatter<Eigen::VectorXi, ConfigurationInfoData>
    VectorXiConfigurationInfoFormatter;

typedef Generic2DDatumFormatter<Eigen::MatrixXd, ConfigurationInfoData>
    MatrixXdConfigurationInfoFormatter;

typedef Generic2DDatumFormatter<Eigen::MatrixXi, ConfigurationInfoData>
    MatrixXiConfigurationInfoFormatter;

ConfigurationInfoFormatter<std::string> configuration_name() {
  return ConfigurationInfoFormatter<std::string>(
      "configuration_name",
      "A name given to the configuration. Consists of the supercell name and "
      "an index. Requires a CASM project.",
      [](ConfigurationInfoData const &data) -> std::string {
        return data.configuration.name();
      });
}

template <typename VectorType>
void write_vector(jsonParser &j, std::vector<std::string> name,
                  VectorType const &v) {
  for (Index i = 0; i < v.size(); ++i) {
    std::string key = name[i];
    if (!j.contains(key)) {
      j[key].put_array();
    }
    j[key].push_back(v(i));
  }
}

template <typename VectorType>
void write_vector(jsonParser &j, std::string name, VectorType const &v) {
  for (Index i = 0; i < v.size(); ++i) {
    std::string key = name + std::to_string(i);
    if (!j.contains(key)) {
      j[key].put_array();
    }
    j[key].push_back(v(i));
  }
}

ConfigurationInfoFormatter<jsonParser> slice() {
  return ConfigurationInfoFormatter<jsonParser>(
      "slice",
      "Return sites within a particular distance of a plane. Requires...",
      [](ConfigurationInfoData const &data) -> jsonParser {
        ParentInputParser parser{data.params_json};
        Eigen::Vector3d b1;
        parser.require(b1, "b1");
        Eigen::Vector3d b2;
        parser.require(b2, "b2");
        Eigen::Vector3d point;
        parser.require(point, "point");
        double distance = 0.1;
        parser.optional(distance, "distance");
        Eigen::Matrix3l T = Eigen::Matrix3l::Identity();
        parser.optional(T, "transformation_matrix_to_supercell");
        std::runtime_error error_if_invalid{"Error reading \"slice\" input"};
        report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

        b1.normalize();
        b2 = b2 - b2.dot(b1) * b1;
        b2.normalize();
        Eigen::Vector3d b3 = b1.cross(b2);
        Eigen::Matrix3d basis = Eigen::Matrix3d::Identity();
        basis.col(0) = b1;
        basis.col(1) = b2;
        basis.col(2) = b3;

        Eigen::Matrix3d basis_inv = basis.transpose()
                                        .colPivHouseholderQr()
                                        .solve(Eigen::Matrix3d::Identity())
                                        .transpose();

        jsonParser j;

        // Make a super-configuration
        Supercell const &supercell = data.configuration.supercell();
        auto shared_supercell = std::make_shared<Supercell const>(
            supercell.shared_prim(),
            supercell.sym_info().transformation_matrix_to_super() * T);

        Configuration superconfig =
            fill_supercell(data.configuration, shared_supercell);

        // Write info about a slice of the super-config
        ConfigDoF const &configdof = superconfig.configdof();
        for (Index l = 0; l < shared_supercell->num_sites(); ++l) {
          Eigen::VectorXd r = shared_supercell->coord(l).const_cart();
          Eigen::VectorXd r_proj = basis_inv * (r - point);
          if (std::abs(r_proj(2)) < distance) {
            // include point & point under change of basis
            write_vector(j, {"rx", "ry", "rz"}, r);
            write_vector(j, {"r1", "r2", "r3"}, r_proj);
            UnitCellCoord bijk = shared_supercell->uccoord(l);
            if (!j.contains("b")) {
              j["b"].put_array();
            }
            j["b"].push_back(bijk.sublattice());
            write_vector(j, {"i", "j", "k"}, bijk.unitcell());

            if (configdof.occupation().size()) {
              if (!j.contains("occ")) {
                j["occ"].put_array();
              }
              j["occ"].push_back(superconfig.occ(l));
              if (!j.contains("species")) {
                j["species"].put_array();
              }
              j["species"].push_back(superconfig.mol(l).name());
            }
            if (!configdof.local_dofs().empty()) {
              for (auto const &local_dof : configdof.local_dofs()) {
                write_vector(j, local_dof.first,
                             local_dof.second.standard_values().col(l));
              }
            }
          }
        }

        if (!configdof.global_dofs().empty()) {
          for (auto const &global_dof : configdof.global_dofs()) {
            write_vector(j, global_dof.first,
                         global_dof.second.standard_values());
          }
        }

        return j;
      });
}

}  // namespace

namespace adapter {

template <typename ToType, typename FromType>
struct Adapter;

template <>
struct Adapter<Configuration, ConfigurationInfoData> {
  Configuration const &operator()(
      ConfigurationInfoData const &adaptable) const {
    return adaptable.configuration;
  }
};

}  // namespace adapter

namespace {

DataFormatterDictionary<ConfigurationInfoData> make_configuration_info_dict() {
  DataFormatterDictionary<ConfigurationInfoData> configuration_info_dict;

  // properties that require ConfigurationInfoData
  configuration_info_dict.insert(configuration_name(), slice());

  // properties that only require configuration
  auto configuration_dict = make_dictionary<Configuration>();
  for (auto it = configuration_dict.begin(); it != configuration_dict.end();
       ++it) {
    if (it->type() == DatumFormatterClass::Property) {
      configuration_info_dict.insert(
          make_datum_formatter_adapter<ConfigurationInfoData, Configuration>(
              *it));
    }
  }

  return configuration_info_dict;
}

}  // namespace

std::string ConfigurationInfoInterface::desc() const {
  std::string description =
      "Get information about a configuration.                             \n\n";

  std::string custom_options =
      "  prim: JSON object (optional, default=prim of current project)    \n"
      "    See `casm format --prim` for details on the prim format.       \n\n"

      "  configuration: configuration                                     \n"
      "    Configuration to get information from.                         \n\n"

      "  params: JSON object (optional, default={})                       \n"
      "    JSON object containing additional parameters for methods that  \n"
      "    require them. See method descriptions for details of how they  \n"
      "    are used or other specific parameters not documented here.     \n\n"

      "  properties: array of string                                      \n"
      "    An array of strings specifying which configuration properties  \n"
      "    to output. The allowed options are:                            \n\n";

  std::stringstream ss;
  auto dict = make_configuration_info_dict();
  dict.print_help(ss, DatumFormatterClass::Property);

  return name() + ": \n\n" + description + custom_options + ss.str();
}

std::string ConfigurationInfoInterface::name() const {
  return "ConfigurationInfo";
}

/// Run `configuration` info method
void ConfigurationInfoInterface::run(jsonParser const &json_options,
                                     PrimClex const *primclex,
                                     fs::path root) const {
  Log &log = CASM::log();

  ParentInputParser parser{json_options};
  std::runtime_error error_if_invalid{"Error reading ConfigurationInfo input"};

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
    if (!root.empty()) {
      try {
        shared_prim = open_shared_prim(root);
      } catch (std::exception &e) {
        parser.insert_error("prim", e.what());
      }
    } else {
      std::stringstream msg;
      msg << "Error in ConfigurationInfo: No \"prim\" in input and no project "
             "provided or found.";
      parser.insert_error("prim", msg.str());
    }
  }
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  // read "configuration"
  std::unique_ptr<Configuration> configuration_ptr =
      parser.require<Configuration>("configuration", shared_prim);

  // read "properties"
  std::vector<std::string> properties;
  parser.require(properties, "properties");

  // read "params"
  jsonParser params_json = jsonParser::object();
  parser.optional(params_json, "params");
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  Configuration const &configuration = *configuration_ptr;
  DataFormatterDictionary<ConfigurationInfoData> configuration_info_dict =
      make_configuration_info_dict();

  // format
  auto formatter = configuration_info_dict.parse(properties);
  jsonParser json;
  ConfigurationInfoData data{primclex, root, shared_prim, configuration,
                             params_json};
  Supercell const &supercell = configuration.supercell();
  if (!supercell.has_primclex()) {
    supercell.set_primclex(&data.primclex());
  }
  formatter.to_json(data, json);
  log << json << std::endl;
}

}  // namespace CASM
