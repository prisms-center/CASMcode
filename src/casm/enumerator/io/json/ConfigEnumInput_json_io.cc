#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/io/json/ConfigDoF_json_io.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "casm/clusterography/io/json/ClusterSpecs_json_io.hh"
#include "casm/crystallography/io/SuperlatticeEnumeratorIO.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"
#include "casm/database/ScelDatabaseTools.hh"
#include "casm/database/io/json_io_impl.hh"
#include "casm/enumerator/ClusterSitesSelector_impl.hh"
#include "casm/enumerator/ConfigEnumInput_impl.hh"
#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"

namespace CASM {

  /// Output ConfigEnumInput to JSON
  jsonParser &to_json(ConfigEnumInput const &config_enum_input, jsonParser &json) {
    json["supercell"] = config_enum_input.configuration().supercell().transf_mat();
    json["configdof"] = config_enum_input.configuration().configdof();
    json["sites"] = config_enum_input.sites();
    return json;
  }

  /// Read ConfigEnumInput from JSON
  ConfigEnumInput jsonConstructor<ConfigEnumInput>::from_json(
    const jsonParser &json,
    std::shared_ptr<Structure const> const &shared_prim,
    DB::Database<Supercell> &supercell_db) {

    jsonParser tjson {json};
    InputParser<ConfigEnumInput> parser {tjson, shared_prim, supercell_db};

    std::runtime_error error_if_invalid {"Error reading ConfigEnumInput from JSON"};
    report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

    return *parser.value;
  }

  /// Read ConfigEnumInput from JSON
  void parse(
    InputParser<ConfigEnumInput> &parser,
    std::shared_ptr<Structure const> const &shared_prim,
    DB::Database<Supercell> &supercell_db) {

    Eigen::Matrix3l T;
    parser.require(T, "supercell");
    auto configdof_ptr = parser.require<ConfigDoF>("configdof", *shared_prim);
    auto sites_ptr = parser.require<std::set<Index>>("site");

    if(parser.valid()) {
      auto supercell_it = make_canonical_and_insert(shared_prim, T, supercell_db).first;
      Configuration configuration {*supercell_it, jsonParser {}, *configdof_ptr};
      parser.value = notstd::make_unique<ConfigEnumInput>(configuration, *sites_ptr);
    }
  }

  /// Read std::vector<ConfigEnumInput> from JSON input, allowing queries from databases
  ///
  /// Note: See `parse` for JSON documentation
  void from_json(
    std::vector<std::pair<std::string, ConfigEnumInput>> &config_enum_input,
    jsonParser const &json,
    std::shared_ptr<Structure const> shared_prim,
    PrimClex const *primclex,
    DB::Database<Supercell> &supercell_db,
    DB::Database<Configuration> &configuration_db) {

    jsonParser tjson {json};
    InputParser<std::vector<std::pair<std::string, ConfigEnumInput>>> parser {tjson, shared_prim, primclex, supercell_db, configuration_db};

    std::runtime_error error_if_invalid {"Error reading std::vector<ConfigEnumInput> from JSON"};
    report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
  }

  /// Make a ScelEnumProps object from JSON input
  void parse(
    InputParser<xtal::ScelEnumProps> &parser,
    DB::Database<Supercell> &supercell_db) {
    int min, max;
    std::string dirs;
    std::string default_dirs {"abc"};
    Eigen::Matrix3i generating_matrix;
    Eigen::Matrix3i default_matrix {Eigen::Matrix3i::Identity()};
    parser.require(min, "min");
    parser.require(max, "max");
    parser.optional_else(dirs, "dirs", default_dirs);

    // support "unit_cell" == supercell name
    if(parser.self.contains("unit_cell") && parser.self["unit_cell"].is_string()) {
      std::string scelname = parser.self["unit_cell"].get<std::string>();
      auto supercell_it = supercell_db.find(scelname);
      if(supercell_it == supercell_db.end()) {
        std::stringstream msg;
        msg << "Error parsing \"unit_cell\": string value not equal to any existing supercell name";
        parser.error.insert(msg.str());
      }
      else {
        generating_matrix = supercell_it->sym_info().transformation_matrix_to_super().cast<int>();
      }
    }
    else {
      parser.optional_else(generating_matrix, "unit_cell", default_matrix);
    }

    if(parser.self.contains("existing_only")) {
      std::stringstream msg;
      msg << "The option \"existing_only\" is no longer supported. "
          << "To select all existing supercells use \"supercell_selection\": \"ALL\".";
      parser.error.insert(msg.str());
    }

    parser.value = notstd::make_unique<xtal::ScelEnumProps>(min, max + 1, dirs, generating_matrix);
  }

  /// Parse JSON to construct initial states for enumeration (as std::vector<std::pair<std::string, ConfigEnumInput>>)
  ///
  /// The result pairs a name, describing the Supercell or Configuration and selected sites name, to
  /// a ConfigEnumInput.
  ///
  /// This method enables several options for specifying initial states for enumeration. There are
  /// two main categories of options: i) options specifying supercells and configurations, and ii)
  /// options restricting which sites the enumeration takes place on.
  ///
  /// Specifying supercells and configurations: These options are all additive in the sense that
  /// the result will be the combination of all that are found. Selecting a Supercell as the initial
  /// state of enumeration is equivalent to selecting the configuration in that supercell with all
  /// values of DoF set to 0.
  ///
  /// Options for specifying supercells and configurations are:
  ///     supercell_selection: string (optional)
  ///         Name of a selection of supercells to use as initial states for enumeration
  ///     scelnames: array of string (optional)
  ///         Array of names of supercell to use as initial states for enumeration
  ///     config_selection: string (optional)
  ///         Name of a selection of configurations to use as initial states for enumeration
  ///     confignames: array of string (optional)
  ///         Array of names of configuration to use as initial states for enumeration
  ///     supercells: object (optional)
  ///         Specifies parameters for enumerating supercells. Options are:
  ///
  ///         min: int (optional, default=1)
  ///             Minimum volume supercells to enumerate
  ///         max: int (required)
  ///             Maximum volume supercells to enumerate
  ///         dirs: string (optional, default="abc")
  ///             Which lattice vectors of unit cell to enumerate over
  ///         unit_cell: 3x3 matrix of int (optional, default=identity matrix)
  ///             The unit cell to tile into supercells.
  ///
  /// Specifying all sites, particular sublattices, particular sites, or particular clusters of
  /// sites to enumerate local DoF. The default behavior if none of these options are
  /// given is selecting all sites for enumeration. The values of DoF on all the sites that are not
  /// selected are frozen.
  ///
  /// Options for selecting sites are:
  ///     sublats: array of integer (optional)
  ///         Indices of sublattices to allow enumeration on. May be used with "sites".
  ///     sites: array of array of integer (optional)
  ///         Indices of sites to allow enumeration on, using [b, i, j, k] notation (b=sublattice
  ///         index, (i,j,k)=unit cell indices). May be used with "sublats". Example:
  ///
  ///             "sites": [
  ///               [0, 0, 0, 0],
  ///               [0, 1, 0, 0],
  ///               [1, 0, 0, 0]
  ///             ]
  ///
  ///     cluster_specs: object (optional)
  ///         JSON object specifying orbits of clusters to generate. Each orbit prototype is used
  ///         to select sites to enumerate on each selected supercell or configuration. If there are
  ///         4 supercells or configurations selected, and there are 10 orbits generated, then there
  ///         will be 4*10=40 ConfigEnumInput generated. The "cluster_specs" option cannot be used
  ///         with the "sublats" or "sites" options.
  ///
  ///
  /// Note:
  /// - Parameter `PrimClex const *primclex` is used while transitioning Supercell to no longer need
  ///   a `PrimClex const *`
  /// - An error is inserted if no ConfigEnumInput are generated
  void parse(
    InputParser<std::vector<std::pair<std::string, ConfigEnumInput>>> &parser,
    std::shared_ptr<Structure const> shared_prim,
    PrimClex const *primclex,
    DB::Database<Supercell> &supercell_db,
    DB::Database<Configuration> &configuration_db) {

    // check for "config_selection" and "confignames"
    DB::Selection<Configuration> config_selection;
    try {
      config_selection = DB::make_selection<Configuration>(
                           configuration_db, parser.self, "confignames", "config_selection");
    }
    catch(std::exception &e) {
      std::stringstream msg;
      msg << "Error creating enumerator initial states from configurations: " << e.what();
      parser.error.insert(msg.str());
    }

    // check for "supercell_selection" and "scelnames"
    DB::Selection<Supercell> supercell_selection;
    try {
      supercell_selection = DB::make_selection<Supercell>(
                              supercell_db, parser.self, "scelnames", "supercell_selection");
    }
    catch(std::exception &e) {
      std::stringstream msg;
      msg << "Error creating enumerator initial states from supercells: " << e.what();
      parser.error.insert(msg.str());
    }

    // check for "supercells"
    auto scel_enum_props_subparser = parser.subparse_if<xtal::ScelEnumProps>("supercells");

    // check for "sublats"
    std::vector<Index> sublats;
    parser.optional(sublats, "sublats");
    for(Index b : sublats) {
      if(b < 0 || b >= shared_prim->basis().size()) {
        std::stringstream msg;
        msg << "Error reading sublats: value out of range [0, " << shared_prim->basis().size() << ")";
        parser.error.insert(msg.str());
      }
    }

    // check for "sites"
    std::vector<UnitCellCoord> sites;
    parser.optional(sites, "sites");
    Index i = 0;
    for(UnitCellCoord site_uccoord : sites) {
      Index b = site_uccoord.sublattice();
      if(b < 0 || b >= shared_prim->basis().size()) {
        std::stringstream msg;
        msg << "Error reading sites[" << i << "]: sublattice index out of range [0, "
            << shared_prim->basis().size() << ")";
        parser.error.insert(msg.str());
      }
      ++i;
    }

    // Do not allow "cluster_specs" with "sublats" or "sites"
    //  (when would this be useful? which order to apply?)
    if((sublats.size() || sites.size()) && parser.self.contains("cluster_specs")) {
      std::stringstream msg;
      msg << "Error creating enumerator initial states: "
          << "cannot include \"cluster_specs\" with \"sublats\" or \"sites\"";
      parser.error.insert(msg.str());
    }

    // check for "cluster_specs"
    auto cluster_specs_subparser = parser.subparse_if<ClusterSpecs>(
                                     "cluster_specs",
                                     shared_prim,
                                     shared_prim->factor_group());
    if(cluster_specs_subparser->value) {
      if(cluster_specs_subparser->value->periodicity_type() != CLUSTER_PERIODICITY_TYPE::PRIM_PERIODIC) {
        std::stringstream msg;
        msg << "Error creating enumerator initial states: "
            << "\"cluster_specs\" method must be \"periodic_max_length\"";
        cluster_specs_subparser->error.insert(msg.str());
      }
    }

    // at this point we have parsed everything except "cluster_specs",
    //   which is parsed later for each ConfigEnumInput
    if(!parser.valid()) {
      return;
    }

    // Use the supercells and configurations input to construct ConfigEnumInput,
    //   by default these have all sites selected
    std::vector<std::pair<std::string, ConfigEnumInput>> config_enum_input;
    for(const auto &config : config_selection.selected()) {
      config_enum_input.emplace_back(config.name(), config);
    }
    for(const auto &scel : supercell_selection.selected()) {
      config_enum_input.emplace_back(scel.name(), scel);
    }
    if(scel_enum_props_subparser->value) {
      ScelEnumByProps enumerator {shared_prim, *scel_enum_props_subparser->value};
      for(auto const &supercell : enumerator) {
        auto supercell_it = supercell_db.insert(supercell).first;
        config_enum_input.emplace_back(supercell_it->name(), *supercell_it);
      }
    }

    // select sublattices and individual sites
    if(sublats.size() || sites.size()) {
      for(auto &input_name_value_pair : config_enum_input) {

        // select sublats & sites, and modify name with sublats & sites
        input_name_value_pair.second.clear_sites();

        std::stringstream ss;
        if(sublats.size()) {
          ss << ", sublats=" << jsonParser {sublats};
          input_name_value_pair.second.select_sublattices(sublats);
        }
        if(sites.size()) {
          ss << ", sites=" << jsonParser {sites};
          input_name_value_pair.second.select_sites(sites);
        }
        input_name_value_pair.first += ss.str();
      }
    }

    // select clusters
    if(cluster_specs_subparser->value != nullptr) {

      // generate orbits from cluster_specs
      auto orbits = cluster_specs_subparser->value->make_periodic_orbits(log());

      // this will generate more ConfigEnumInput, with cluster sites selected
      std::vector<std::pair<std::string, ConfigEnumInput>> all_with_cluster_sites;
      for(auto &input_name_value_pair : config_enum_input) {
        std::vector<ConfigEnumInput> with_cluster_sites;
        select_cluster_sites(input_name_value_pair.second, orbits, std::back_inserter(with_cluster_sites));

        // modify name with cluster site indices
        for(auto const &value : with_cluster_sites) {
          std::stringstream ss;
          ss << ", cluster_sites=" << jsonParser {value.sites()};
          std::string name = input_name_value_pair.first + ss.str();
          all_with_cluster_sites.emplace_back(name, value);
        }
      }
      config_enum_input = std::move(all_with_cluster_sites);
    }

    // Move all constructed ConfigEnumInput into the parser.value
    parser.value = notstd::make_unique<std::vector<std::pair<std::string, ConfigEnumInput>>>(std::move(config_enum_input));

    if(parser.value->size() == 0) {
      std::stringstream msg;
      msg << "Error creating enumerator initial states: No supercells or configurations.";
      parser.error.insert(msg.str());
    }
  }

}
