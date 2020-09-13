#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/ScelEnum_impl.hh"
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
    std::vector<ConfigEnumInput> &config_enum_input,
    jsonParser const &json,
    std::shared_ptr<Structure const> shared_prim,
    DB::Database<Supercell> &supercell_db,
    DB::Database<Configuration> &configuration_db) {

    jsonParser tjson {json};
    InputParser<std::vector<ConfigEnumInput>> parser {tjson, shared_prim, supercell_db, configuration_db};

    std::runtime_error error_if_invalid {"Error reading std::vector<ConfigEnumInput> from JSON"};
    report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
  }


  /// Parse JSON to construct initial states for enumeration (as std::vector<ConfigEnumInput>)
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
  /// Options are:
  ///     supercell_selection: string (optional)
  ///         Name of a selection of supercells to use as initial states for enumeration
  ///     scelnames: array of string (optional)
  ///         Array of names of supercell to use as initial states for enumeration
  ///     config_selection: string (optional)
  ///         Name of a selection of configurations to use as initial states for enumeration
  ///     confignames: array of string (optional)
  ///         Array of names of configuration to use as initial states for enumeration. Also accepts
  ///         "scelnames" for this option.
  ///     supercells: object (optional)
  ///         Specifies parameters for enumerating supercells. Options are:
  ///
  ///         min: int (required)
  ///             Minimum volume supercells to enumerate
  ///         max: int (required)
  ///             Maximum volume supercells to enumerate
  ///         dirs: string (optional, default="abc")
  ///             Which lattice vectors of unit cell to enumerate over
  ///         unit_cell: 3x3 matrix of int (optional, default=identity matrix)
  ///             The unit cell to tile into supercells.
  ///
  /// Specifying all sites, particular sublattices, particular sites, or particular clusters of
  /// sites to enumerate local DoF: These options are all additive in the sense that the result will
  /// be the combination of all that are found. The default behavior if none of these options are
  /// given is selecting all sites for enumeration. The values of DoF on all the sites that are not
  /// selected are frozen.
  ///
  /// Options are:
  ///     sublats: array of integer (optional)
  ///         Indices of sublattices to allow enumeration on.
  ///     sites: array of array of integer (optional)
  ///         Indices of sites to allow enumeration on, using [b, i, j, k] notation (b=sublattice
  ///         index, (i,j,k)=unit cell indices). Example:
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
  ///         will be 4*10=40 ConfigEnumInput generated.
  ///
  void parse(
    InputParser<std::vector<ConfigEnumInput>> &parser,
    std::shared_ptr<Structure const> shared_prim,
    DB::Database<Supercell> &supercell_db,
    DB::Database<Configuration> &configuration_db) {

    // check for "config_selection" and "confignames"
    DB::Selection<Configuration> config_selection;
    if(parser.self.contains("confignames") || parser.self.contains("config_selection")) {
      try {
        config_selection = DB::make_selection<Configuration>(
                             configuration_db, parser.self, "confignames", "config_selection");
      }
      catch(std::exception &e) {
        parser.error.insert(std::string("Error creating enumerator initial states from configurations: ") + e.what());
      }
    }

    // check for "supercell_selection" and "scelnames"
    DB::Selection<Supercell> supercell_selection;
    if(parser.self.contains("scelnames") || parser.self.contains("supercell_selection")) {
      try {
        supercell_selection = DB::make_selection<Supercell>(
                                supercell_db, parser.self, "scelnames", "supercell_selection");
      }
      catch(std::exception &e) {
        parser.error.insert(std::string("Error creating enumerator initial states from supercells: ") + e.what());
      }
    }

    // check for "supercells"
    auto scel_enum_props_subparser = parser.subparse_if<xtal::ScelEnumProps>("supercells");

    // check for "sublats"
    std::vector<Index> sublats;
    parser.optional(sublats, "sublats");
    std::cout << "sublats: " << sublats << std::endl;
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
    jsonParser tmp;
    tmp["sites"] = sites;
    std::cout << "sites: \n" << tmp << std::endl;
    Index i = 0;
    for(UnitCellCoord site_uccoord : sites) {
      Index b = site_uccoord.sublattice();
      if(b < 0 || b >= shared_prim->basis().size()) {
        std::stringstream msg;
        msg << "Error reading sites[" << i << "]: sublattice index out of range [0, " << shared_prim->basis().size() << ")";
        parser.error.insert(msg.str());
      }
      ++i;
    }

    // check for "cluster_specs"
    std::cout << "hi 0" << std::endl;
    auto cluster_specs_subparser = parser.subparse_if<ClusterSpecs>("cluster_specs", shared_prim);

    std::cout << "hi 1" << std::endl;
    if(!parser.valid()) {
      return;
    }

    // Use the parsed input to construct ConfigEnumInput:
    std::cout << "hi 2" << std::endl;
    parser.value = notstd::make_unique<std::vector<ConfigEnumInput>>();
    auto &config_enum_input = *parser.value;

    std::cout << "hi 3" << std::endl;
    for(const auto &config : config_selection.selected()) {
      config_enum_input.emplace_back(config);
    }
    std::cout << "hi 4" << std::endl;
    for(const auto &scel : supercell_selection.selected()) {
      config_enum_input.emplace_back(scel);
    }
    std::cout << "hi 5" << std::endl;
    if(scel_enum_props_subparser->value) {
      ScelEnumByProps enumerator {shared_prim, *scel_enum_props_subparser->value};
      for(const auto &scel : enumerator) {
        auto scel_it = supercell_db.insert(scel).first;
        config_enum_input.emplace_back(*scel_it);
      }
    }
    std::cout << "hi 6" << std::endl;
    if(sublats.size() || sites.size()) {
      for(ConfigEnumInput &input : config_enum_input) {
        input.clear_sites();
        input.select_sublattices(sublats);
        input.select_sites(sites);
      }
    }
    std::cout << "hi 7" << std::endl;
    if(cluster_specs_subparser->value) {
      ClusterSpecs const &cluster_specs = *cluster_specs_subparser->value;
      std::vector<ConfigEnumInput> with_cluster_sites;
      for(ConfigEnumInput &input : config_enum_input) {
        auto selector = make_cluster_sites_selector(input, std::back_inserter(with_cluster_sites));
        for_all_orbits(cluster_specs, null_log(), selector);
      }
      config_enum_input = std::move(with_cluster_sites);
    }
    std::cout << "hi 8" << std::endl;

  }

}
