#include "casm/app/info/methods/ConfigurationInfoInterface.hh"

#include "casm/app/ProjectSettings.hh"
#include "casm/app/info/InfoInterface_impl.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/casm_io/dataformatter/DataFormatterTools_impl.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/dataformatter/DatumFormatterAdapter.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/optional.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/FillSupercell.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/io/json/Configuration_json_io.hh"
#include "casm/clex/io/stream/Configuration_stream_io.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/database/ConfigDatabaseTools_impl.hh"
#include "casm/database/DatabaseTypes_impl.hh"
#include "casm/enumerator/DoFSpace.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"

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
                        std::optional<Configuration> _optional_configuration,
                        jsonParser const &_params_json)
      : shared_prim(_shared_prim),
        optional_configuration(_optional_configuration),
        params_json(_params_json),
        m_primclex(_primclex),
        m_root(_root) {}

  std::shared_ptr<Structure const> shared_prim;

  std::optional<Configuration> optional_configuration;

  jsonParser const &params_json;

  Configuration const &configuration() const {
    if (!optional_configuration.has_value()) {
      throw std::runtime_error("Error in ConfigurationInfo: no configuration");
    }
    return *optional_configuration;
  }

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
        return data.configuration().name();
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
      "Return sites within a particular distance of a plane in a supercell of "
      "a configuration. Requires the following params: \"b1\" and \"b2\", "
      "(array, size=3) Cartesian vectors in the plane; \"point\", (array, "
      "size=3) a Cartesian coordinate point in the plane; \"distance\", "
      "(number) the distance from the plane to include sites; "
      "\"transformation_matrix_to_supercell\", (2d array, shape=(3,3)), the "
      "transformation matrix for a supercell which is filled by the "
      "configuration. The output is several arrays where the i-th element "
      "describes the i-site within the specified distance to the specified "
      "plane. Information includes: \"rx\", \"ry\", \"rz\", Cartesian "
      "coordinates; \"r1\", \"r2\", \"r3\", projected coordinates, where "
      "\"r1\" and \"r2\" are coordinates lieing along normalized \"b1\" and "
      "\"b2\" axes, and \"r3\" is the distance from the plane; \"b\", \"i\", "
      "\"j\", \"k\", are the integral coordinates of the site; \"occ\", the "
      "indicator value of the sites occupant; \"species\", the name of the "
      "occupant species; \"<local_dof_name><j>\", the j-th component of local "
      "DoF on the site, in the standard basis; \"<global_dof_name><j>\", the "
      "j-th component of global DoF (constant for all sites).",
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
        Supercell const &supercell = data.configuration().supercell();
        auto shared_supercell = std::make_shared<Supercell const>(
            supercell.shared_prim(),
            supercell.sym_info().transformation_matrix_to_super() * T);

        Configuration superconfig =
            fill_supercell(data.configuration(), shared_supercell);

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

VectorXdConfigurationInfoFormatter normal_coordinate() {
  return VectorXdConfigurationInfoFormatter(
      "normal_coordinate",
      "Return DoF values expressed as a normal coordinate in a DoF space "
      "basis, for instance as generated by the `--dof-space-analysis` or "
      "`--config-space-analysis` methods of `casm sym`. This may be used for "
      "global DoF with any configuration. This may be used for occ or local "
      "DoF for configurations which have the same supercell as the DoF space. "
      "Otherwise, use `normal_coordinate_at` or `order_parameter`. For local "
      "DoF, the DoF space may be a subset of sites. Requires the following "
      "params: \"dof_space\", (object) a DoF space basis description, as "
      "output by `--dof-space-analysis` or `--config-space-analysis` methods, "
      "including: \"dof\" (string) the DoF type; \"basis\" (2d array, rows are "
      "axes and must have dimension matching the DoF space dimension, the "
      "number of columns must be less than or equal to the DoF space "
      "dimension), the DoF space basis (default is identity); "
      "\"transformation_matrix_to_supercell\" (2d array, shape=(3,3)) the "
      "matrix specifying the supercell lattice vectors (for occ and local "
      "DoF); and \"sites\", (array of int), specifying the linear indices of "
      "sites in the supercell whose DoF are inclued in the DoF space (for occ "
      "and local DoF). The output is a vector, of size equal to the DoF space "
      "basis size.",
      [](ConfigurationInfoData const &data) -> Eigen::VectorXd {
        ParentInputParser parser{data.params_json};
        std::unique_ptr<DoFSpace> dof_space_ptr =
            parser.require<DoFSpace>("dof_space", data.shared_prim);
        std::runtime_error error_if_invalid{
            "Error reading \"normal_coordinate\" input"};
        report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

        DoFSpace const &dof_space = *dof_space_ptr;
        return get_normal_coordinate(data.configuration(), dof_space);
      });
}

ConfigurationInfoFormatter<jsonParser> normal_coordinate_at_each_unitcell() {
  return ConfigurationInfoFormatter<jsonParser>(
      "normal_coordinate_at_each_unitcell",
      "Return DoF values expressed as a normal coordinate in a DoF space "
      "basis, calculated for each unitcell in a configuration. This method is "
      "similar to \"normal_coordinate\", but for occ DoF and local DoF only, "
      "and it acts on a subset of configuration, so the configuration does not "
      "need to have the same supercell as the DoF space. Requires the "
      "following params: \"dof_space\" (see \"normal_coordinate\"); "
      "\"unitcell\" (array of int, size=3) the (i,j,k) coordinates of the unit "
      "cell at which the DoF are read from the configuration. The output is "
      "two arrays: \"unitcell\" an array of (i,j,k) indices for each unit cell "
      "in the configuration; \"normal_coordinate\" an array of normal "
      "coordinates calculated at the corresponding unit cell.",
      [](ConfigurationInfoData const &data) -> jsonParser {
        ParentInputParser parser{data.params_json};
        std::unique_ptr<DoFSpace> dof_space_ptr =
            parser.require<DoFSpace>("dof_space", data.shared_prim);
        std::runtime_error error_if_invalid{
            "Error reading \"normal_coordinate_at_each_unitcell\" input"};
        report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

        DoFSpace const &dof_space = *dof_space_ptr;
        DoFSpaceIndexConverter index_converter{data.configuration(), dof_space};
        auto const &unitcell_index_converter = data.configuration()
                                                   .supercell()
                                                   .sym_info()
                                                   .unitcell_index_converter();

        jsonParser json;
        json["unitcell"].put_array();
        json["normal_coordinate"].put_array();
        for (int i = 0; i < unitcell_index_converter.total_sites(); ++i) {
          xtal::UnitCell unitcell = unitcell_index_converter(i);
          Eigen::VectorXd normal_coordinate = get_normal_coordinate_at(
              data.configuration(), dof_space, index_converter, unitcell);
          {
            jsonParser tjson;
            to_json(unitcell, tjson, jsonParser::as_array());
            json["unitcell"].push_back(tjson);
          }
          {
            jsonParser tjson;
            to_json(normal_coordinate, tjson, jsonParser::as_array());
            json["normal_coordinate"].push_back(tjson);
          }
        }
        return json;
      });
}

VectorXdConfigurationInfoFormatter order_parameter() {
  return VectorXdConfigurationInfoFormatter(
      "order_parameter",
      "Return order parameters calculated as the mean value of the normal "
      "coordinate in a DoF space basis, averaged over an entire configuration. "
      "This method is similar to \"normal_coordinate\", but for occ DoF and "
      "local DoF it calculates the average value over the tilings of a "
      "configuration into the commensurate supercell of the configuration and "
      "the DoF space. The resulting value may be useful as an order parameter. "
      "The configuration does not need to have the same supercell as the DoF "
      "space. Requires the following params: \"dof_space\" (see "
      "\"normal_coordinate\"). The output is a vector, of size equal to the "
      "DoF space basis size.",
      [](ConfigurationInfoData const &data) -> Eigen::VectorXd {
        ParentInputParser parser{data.params_json};
        std::unique_ptr<DoFSpace> dof_space_ptr =
            parser.require<DoFSpace>("dof_space", data.shared_prim);
        std::runtime_error error_if_invalid{
            "Error reading \"order_parameter\" input"};
        report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

        DoFSpace const &dof_space = *dof_space_ptr;
        return get_mean_normal_coordinate(data.configuration(), dof_space);
      });
}

ConfigurationInfoFormatter<jsonParser> insert() {
  return ConfigurationInfoFormatter<jsonParser>(
      "insert",
      "Insert configurations into the project configuration database. May take "
      "the following params: \"config_list\" (optional, array of configuration "
      "objects), additional configurations to be imported into the "
      "configuration database; \"enum_output_file\" (optional, path to JSON "
      "file containing an array of with \"config\" attributes), additional "
      "configurations to be imported into the configuration database; "
      "\"insert_primitive_only\" (bool, default=false), if true, make and "
      "insert primitive configurations only; \"dry_run\" (bool) if true, "
      "insert configurations and report results, but do not save database. "
      "The output is an array of JSON objects, giving the configurations "
      "along with the name they were given when inserted.",
      [](ConfigurationInfoData const &data) -> jsonParser {
        auto shared_prim = data.primclex().shared_prim();
        auto &configuration_db = data.primclex().db<Configuration>();
        auto &supercell_db = data.primclex().db<Supercell>();

        ParentInputParser parser{data.params_json};
        std::runtime_error error_if_invalid{"Error reading \"insert\" input"};

        bool insert_primitive_only = false;
        parser.optional(insert_primitive_only, "insert_primitive_only");

        bool dry_run = false;
        parser.optional(dry_run, "dry_run");

        std::set<Configuration> configs;

        if (data.optional_configuration.has_value()) {
          // --- hack to get input configuration with primclex shared_prim in
          // all cases ---
          auto shared_supercell = std::make_shared<Supercell const>(
              shared_prim, data.configuration()
                               .supercell()
                               .sym_info()
                               .transformation_matrix_to_super());
          Configuration configuration(shared_supercell);
          configuration.configdof().values() =
              data.configuration().configdof().values();
          // --- end ---

          configs.insert(configuration);
        }
        if (parser.self.contains("config_list")) {
          int n = parser.self["config_list"].size();
          for (int i = 0; i < n; ++i) {
            std::unique_ptr<Configuration> configuration_ptr =
                parser.require<Configuration>(
                    fs::path{"config_list"} / std::to_string(i), shared_prim);
            report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
            configs.insert(*configuration_ptr);
          }
        }
        if (parser.self.contains("enum_output_file")) {
          std::string enum_output_file;
          parser.require(enum_output_file, "enum_output_file");
          report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

          fs::path filepath{enum_output_file};
          if (!fs::exists(filepath)) {
            parser.insert_error("enum_output_file", "File does not exist");
            report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
          }
          jsonParser json{filepath};
          ParentInputParser enum_output_parser{json};
          std::runtime_error enum_output_error_if_invalid{
              "Error reading \"insert\"/\"enum_output_file\" file"};
          if (!enum_output_parser.self.is_array()) {
            enum_output_parser.error.insert("Not an array");
            report_and_throw_if_invalid(enum_output_parser, CASM::log(),
                                        enum_output_error_if_invalid);
          }
          int n = enum_output_parser.self.size();
          for (int i = 0; i < n; ++i) {
            std::unique_ptr<Configuration> configuration_ptr =
                enum_output_parser.require<Configuration>(
                    fs::path{std::to_string(i)} / "config", shared_prim);
            report_and_throw_if_invalid(enum_output_parser, CASM::log(),
                                        enum_output_error_if_invalid);
            configs.insert(*configuration_ptr);
          }
        }

        std::map<Configuration, std::string> inserted_configs;
        jsonParser json = jsonParser::array();
        for (auto const &config : configs) {
          ConfigInsertResult result = DB::make_canonical_and_insert(
              config, supercell_db, configuration_db, insert_primitive_only);

          jsonParser tjson;
          tjson["config"] = config;
          if (!insert_primitive_only) {
            tjson["is_new"] = result.insert_canonical;
            tjson["canonical_config"] = *result.canonical_it;
            tjson["canonical_name"] = result.canonical_it->name();
          }
          tjson["is_new_primitive"] = result.insert_primitive;
          tjson["primitive_config"] = *result.primitive_it;
          tjson["primitive_name"] = result.primitive_it->name();
          json.push_back(tjson);
        }

        // Set primclex pointer
        for (auto const &config : configuration_db) {
          if (!config.supercell().has_primclex()) {
            config.supercell().set_primclex(&data.primclex());
          }
        }

        if (!dry_run) {
          configuration_db.commit();
          supercell_db.commit();
        }

        return json;
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
    if (!adaptable.optional_configuration.has_value()) {
      throw std::runtime_error("Error in ConfigurationInfo: no configuration");
    }
    return *adaptable.optional_configuration;
  }
};

}  // namespace adapter

namespace {

DataFormatterDictionary<ConfigurationInfoData> make_configuration_info_dict() {
  DataFormatterDictionary<ConfigurationInfoData> configuration_info_dict;

  // properties that require ConfigurationInfoData
  configuration_info_dict.insert(configuration_name(), normal_coordinate(),
                                 normal_coordinate_at_each_unitcell(),
                                 order_parameter(), insert(), slice());

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
  ss << name() + ": \n\n" + description + custom_options;
  auto dict = make_configuration_info_dict();
  print_info_desc(dict, ss);
  return ss.str();
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
  std::optional<Configuration> optional_configuration;
  parser.optional(optional_configuration, "configuration", shared_prim);

  // read "properties"
  std::vector<std::string> properties;
  parser.require(properties, "properties");

  // read "params"
  jsonParser params_json = jsonParser::object();
  parser.optional(params_json, "params");
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  DataFormatterDictionary<ConfigurationInfoData> configuration_info_dict =
      make_configuration_info_dict();

  // format
  auto formatter = configuration_info_dict.parse(properties);
  jsonParser json;
  ConfigurationInfoData data{primclex, root, shared_prim,
                             optional_configuration, params_json};
  if (data.optional_configuration.has_value()) {
    Supercell const &supercell = data.configuration().supercell();
    if (!supercell.has_primclex()) {
      supercell.set_primclex(&data.primclex());
    }
  }
  formatter.to_json(data, json);
  log << json << std::endl;
}

}  // namespace CASM
