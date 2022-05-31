#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/FillSupercell_impl.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/io/json/ConfigDoF_json_io.hh"
#include "casm/clex/io/json/Configuration_json_io.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "casm/clusterography/io/json/ClusterSpecs_json_io.hh"
#include "casm/crystallography/io/SuperlatticeEnumeratorIO.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"
#include "casm/database/ScelDatabaseTools.hh"
#include "casm/database/io/json_io_impl.hh"
#include "casm/enumerator/ClusterSitesSelector_impl.hh"
#include "casm/enumerator/ConfigEnumInput_impl.hh"

namespace CASM {

/// Output ConfigEnumInput to JSON
jsonParser &to_json(ConfigEnumInput const &config_enum_input,
                    jsonParser &json) {
  json["configuration"] = config_enum_input.configuration();
  json["sites"] = config_enum_input.sites();
  return json;
}

/// Read ConfigEnumInput from JSON
ConfigEnumInput jsonConstructor<ConfigEnumInput>::from_json(
    const jsonParser &json,
    std::shared_ptr<Structure const> const &shared_prim) {
  jsonParser tjson{json};
  InputParser<ConfigEnumInput> parser{tjson, shared_prim};

  std::runtime_error error_if_invalid{
      "Error reading ConfigEnumInput from JSON"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

  return *parser.value;
}

/// Read ConfigEnumInput from JSON
void parse(InputParser<ConfigEnumInput> &parser,
           std::shared_ptr<Structure const> const &shared_prim) {
  auto configuration_ptr =
      parser.require<Configuration>("configuration", shared_prim);
  auto sites_ptr = parser.require<std::set<Index>>("sites");

  if (parser.valid()) {
    parser.value =
        notstd::make_unique<ConfigEnumInput>(*configuration_ptr, *sites_ptr);
  }
}

/// Read std::vector<ConfigEnumInput> from JSON input, allowing queries from
/// databases
///
/// Note: See `parse` for JSON documentation
void from_json(
    std::vector<std::pair<std::string, ConfigEnumInput>> &config_enum_input,
    jsonParser const &json, std::shared_ptr<Structure const> shared_prim,
    PrimClex const *primclex, DB::Database<Supercell> &supercell_db,
    DB::Database<Configuration> &configuration_db) {
  jsonParser tjson{json};
  InputParser<std::vector<std::pair<std::string, ConfigEnumInput>>> parser{
      tjson, shared_prim, primclex, supercell_db, configuration_db};

  std::runtime_error error_if_invalid{
      "Error reading std::vector<ConfigEnumInput> from JSON"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
}

/// Make a ScelEnumProps object from JSON input
///
/// Notes:
/// - Default "min" == 1
/// - Default "dirs" == "abc"
/// - Default "unit_cell" == 3x3 identity matrix
void parse(InputParser<xtal::ScelEnumProps> &parser,
           DB::Database<Supercell> &supercell_db) {
  int min, max;
  std::string dirs;
  int default_min = 1;
  std::string default_dirs{"abc"};
  Eigen::Matrix3i generating_matrix;
  Eigen::Matrix3i default_matrix{Eigen::Matrix3i::Identity()};
  bool diagonal_only;
  bool fixed_shape;

  parser.optional_else(min, "min", default_min);
  parser.require(max, "max");
  parser.optional_else(dirs, "dirs", default_dirs);

  // support "unit_cell" == supercell name
  if (parser.self.contains("unit_cell") &&
      parser.self["unit_cell"].is_string()) {
    std::string scelname = parser.self["unit_cell"].get<std::string>();
    auto supercell_it = supercell_db.find(scelname);
    if (supercell_it == supercell_db.end()) {
      std::stringstream msg;
      msg << "Error parsing \"unit_cell\": string value not equal to any "
             "existing supercell name";
      parser.error.insert(msg.str());
    } else {
      generating_matrix =
          supercell_it->sym_info().transformation_matrix_to_super().cast<int>();
    }
  } else {
    parser.optional_else(generating_matrix, "unit_cell", default_matrix);
  }

  parser.optional_else(diagonal_only, "diagonal_only", false);
  parser.optional_else(fixed_shape, "fixed_shape", false);

  if (parser.self.contains("existing_only")) {
    std::stringstream msg;
    msg << "The option \"existing_only\" is no longer supported. "
        << "To select all existing supercells use \"supercell_selection\": "
           "\"ALL\".";
    parser.error.insert(msg.str());
  }

  if (parser.valid()) {
    parser.value = notstd::make_unique<xtal::ScelEnumProps>(
        min, max + 1, dirs, generating_matrix, diagonal_only, fixed_shape);
  }
}

/// A string describing the JSON format for parsing named ConfigEnumInput
std::string parse_ConfigEnumInput_desc() {
  return "  scelnames: array of strings (optional, override with --scelnames)  "
         "               \n"
         "    Names of supercells used as input states. All sites will be set "
         "to the first    \n"
         "    listed occupant, and all other DoFs will be set to zero. Ex:     "
         "               \n\n"

         "        \"scelnames\": [\"SCEL1_1_1_1_0_0_0\", "
         "\"SCEL2_2_1_1_0_0_0\"]              \n\n"

         "  supercell_selection: string (optional)                             "
         "               \n"
         "    Name of a selection of supercells to use as input states.        "
         "               \n\n"

         "  supercells: object, ScelEnum JSON settings (optional, override "
         "with --min, --max) \n"
         "    Indicate supercells to use as input states in terms of size and "
         "unit cell via   \n"
         "    a JSON object conforming to the format of 'ScelEnum' JSON "
         "settings \"min\",     \n"
         "    \"max\", \"dirs\", and \"unit_cell\". See 'ScelEnum' description "
         "for more       \n"
         "    details.                                                         "
         "               \n\n"

         "  confignames: array of strings (optional, override with "
         "--confignames)             \n"
         "    Names of configurations to be used as input states. All "
         "specified sublattices or\n"
         "    sites will be enumerated on and all other DoFs will              "
         "               \n"
         "    maintain the values of the initial state. Ex:                    "
         "               \n\n"

         "        \"confignames\": [\"SCEL1_1_1_1_0_0_0/1\", "
         "\"SCEL2_2_1_1_0_0_0/3\"]        \n\n"

         "  config_selection: string (optional)                                "
         "               \n"
         "    Name of a selection of configurations to use as initial states.  "
         "               \n\n"

         "  config_list: array (optional)                                      "
         "               \n"
         "    Allows direct input of configurations that may not already exist "
         "in the project \n"
         "    database via a JSON array of objects conforming to the "
         "\"config.json\" format.  \n"
         "    If the \"dof\" component is not included all sites will be set "
         "to the first     \n"
         "    listed occupant, and all other DoFs will be set to zero, as if "
         "specifying a     \n"
         "    a supercell. Ex:                                                 "
         "               \n\n"

         "        \"config_list\" : [                                          "
         "               \n"
         "          {                                                          "
         "               \n"
         "            \"transformation_matrix_to_supercell\": [                "
         "               \n"
         "              [ 0, -1, -1 ],                                         "
         "               \n"
         "              [ 0, 1, -1 ],                                          "
         "               \n"
         "              [ 1, 0, 1 ]                                            "
         "               \n"
         "            ],                                                       "
         "               \n"
         "            \"identifier\": \"custom_strain.1\",                     "
         "               \n"
         "            \"dof\": {                                               "
         "               \n"
         "              \"global_dofs\" : {                                    "
         "               \n"
         "                \"Hstrain\" : {                                      "
         "               \n"
         "                  \"values\" : [ 0.100000000000, 0.100000000000, "
         "0.100000000000, 0.000000000000, 0.000000000000, 0.000000000000 ]\n"
         "                }                                                    "
         "               \n"
         "              },                                                     "
         "               \n"
         "              \"occ\" : [ 0, 0 ]                                     "
         "               \n"
         "            }                                                        "
         "               \n"
         "          }                                                          "
         "               \n"
         "        ]                                                            "
         "               \n\n"

         "  fill: JSON object (optional, default none)                       \n"
         "    If present, the \"fill\" option specifies that the input       \n"
         "    configurations specified previously (by the \"confignames\",   \n"
         "    \"config_selection\", \"config_list\", \"scelnames\",          \n"
         "    \"supercell_selection\", and \"supercells\" parameters) should \n"
         "    be used to fill supercells specified by this object. The       \n"
         "    supercells to \"fill\"  may be specified by \"scelnames\",     \n"
         "    \"supercell_selection\", and/or \"supercells\" parameters with \n"
         "    the same format as described above. Any configurations that    \n"
         "    cannot exactly tile the specified supercells are skipped, with \n"
         "    a warning.\n\n"

         "    As an example, in a project with an fcc prim, enumeration of   \n"
         "    perturbations of a configuration in supercells of the          \n"
         "    conventional fcc cell can be done with:                        \n"
         "                                                                   \n"
         "        \"confignames\": [\"SCEL4_2_2_1_1_1_0/1\"],                \n"
         "        \"fill\": {                                                \n"
         "          \"supercells\": {                                        \n"
         "            \"min\": 1,                                            \n"
         "            \"max\": 4,                                            \n"
         "            \"unit_cell\": [                                       \n"
         "              [ -1,  1,  1 ],                                      \n"
         "              [  1, -1,  1 ],                                      \n"
         "              [  1,  1, -1 ]                                       \n"
         "            ]                                                      \n"
         "          }                                                        \n"
         "        },                                                         \n"
         "        \"cluster_specs\": {                                       \n"
         "           ...                                                     \n"
         "        }\n\n"

         "  sublats: array of integers (optional, default none)                "
         "               \n"
         "    Selects sites by specifying sublattices. Each sublattice index "
         "corresponds to a \n"
         "    basis site in prim.json, indexed from 0. Ex:                     "
         "               \n\n"

         "        \"sublats\" : [0, 2]                                         "
         "               \n\n"

         "  sites: array of 4-entry integer arrays (optional, default none) \n"
         "    Selects sites by [b,i,j,k] convention, where 'b' is sublattice "
         "index and [i,j,k]\n"
         "    specifies linear combinations of primitive-cell lattice vectors. "
         "Ex:            \n\n"

         "        \"sites\" : [[0,0,0,0],                                      "
         "               \n"
         "                     [2,0,0,0]]                                      "
         "               \n\n"

         "  cluster_specs: object (optional)                                   "
         "               \n"
         "    JSON object specifying orbits of clusters to generate. Each "
         "orbit prototype is  \n"
         "    used to select sites on each input supercell or configuration. "
         "If there are 4   \n"
         "    supercells or configurations selected, and there are 10 orbits "
         "generated, then  \n"
         "    there will be 4*10=40 initial states generated. The "
         "\"cluster_specs\" option    \n"
         "    cannot be used with the \"sublats\" or \"sites\" options. Expect "
         "format is:     \n\n"

         "    method: string (required)                                        "
         "               \n"
         "      Specify which cluster orbit generating method will be used. "
         "Supported:        \n\n"

         "      - \"periodic_max_length\": Clusters differing by a lattice "
         "translation are    \n"
         "        considered equivalent. Cluster generation is truncated by "
         "specifying the    \n"
         "        maximum distance between sites in a cluster for 2-point, "
         "3-point, etc.      \n"
         "        clusters. The point clusters comprising the asymmetric unit "
         "of the prim     \n"
         "        structure are always included. After the cluster orbits are "
         "generated using \n"
         "        the prim factor group symmetry, the orbits are broken into "
         "sub-orbits       \n"
         "        reflecting the configuration factor group symmetry of the "
         "input state. Any  \n"
         "        orbits that are duplicated under periodic boundary "
         "conditions are removed.  \n\n"

         "    params: object (required) \n"
         "      Specifies parameters for the method selected by `method`. "
         "Options depend on the \n"
         "      `method` chosen: \n\n"

         "      For method==\"periodic_max_length\": \n"
         "        orbit_branch_specs: object (optional)\n"
         "          Cluster generation is truncated by specifying the maximum "
         "distance\n"
         "          between sites in a cluster for each orbit branch (i.e. "
         "2-point, 3-point, etc.\n"
         "          clusters). The 1-point clusters comprising the asymmetric "
         "unit of the prim\n"
         "          structure are always included. \n\n"

         "            Example: \n"
         "              \"orbit_branch_specs\": {\n"
         "                \"2\": { \"max_length\": 10.0 },\n"
         "                \"3\": { \"max_length\": 8.0 },\n"
         "                       ...\n"
         "              }\n\n"

         "        orbit_specs: array (optional) \n"
         "          An array of clusters which are used to generate and "
         "include orbits of clusters \n"
         "          whether or not they meet the `max_length` truncation "
         "criteria. See the \n"
         "          cluster input format below. Use the "
         "\"include_subclusters\" option to force \n"
         "          generation of orbits for all subclusters of the specified "
         "cluster. \n"

         "            Example cluster, with \"Direct\" coordinates: \n"
         "              { \n"
         "                \"coordinate_mode\" : \"Direct\", \n"
         "                \"sites\" : [ \n"
         "                  [ 0.000000000000, 0.000000000000, 0.000000000000 "
         "], \n"
         "                  [ 1.000000000000, 0.000000000000, 0.000000000000 "
         "], \n"
         "                  [ 2.000000000000, 0.000000000000, 0.000000000000 "
         "], \n"
         "                  [ 3.000000000000, 0.000000000000, 0.000000000000 "
         "]], \n"
         "                \"include_subclusters\" : true \n"
         "              } \n\n"

         "            Example cluster, with \"Integral\" coordinates: \n"
         "              { \n"
         "                \"coordinate_mode\" : \"Integral\", \n"
         "                \"sites\" : [ \n"
         "                  [ 0, 0, 0, 0 ], \n"
         "                  [ 0, 1, 0, 0 ], \n"
         "                  [ 1, 0, 0, 0 ]], \n"
         "                \"include_subclusters\" : true \n"
         "              } \n\n"

         "      For method==\"local_max_length\": \n"
         "        phenomenal: object (required)\n"
         "          The \"phenomenal\" cluster about which local clusters are\n"
         "          generated. See the cluster input format below.\n\n"

         "        generating_group: array of int (required)\n"
         "          An array of symop indices into the prim structure factor \n"
         "          group specifying the invariant group of the              \n"
         "          \"phenomenal\" cluster which should be used to for       \n"
         "          generating local cluster orbits. In some contexts, the   \n"
         "          relevant symmetry is lower than that determined from the \n"
         "          phenomenal cluster sites alone.\n\n"

         "        orbit_branch_specs: object (required)\n"
         "          All sites within `cutoff_radius` distance of any site in \n"
         "          the phenomenal cluster are considered candidates for     \n"
         "          inclusion in clusters of a particular orbit branch.      \n"
         "          Cluster generation is truncated by specifying the maximum\n"
         "          distance between sites in a cluster for each orbit       \n"
         "          branch. The `max_length` parameter is not necessary for  \n"
         "          1-point clusters and ignored if present.\n\n"

         "            Example:                                               \n"
         "              \"orbit_branch_specs\": {                            \n"
         "                \"1\": { \"cutoff_radius\": 6.0 },                 \n"
         "                \"2\": {                                           \n"
         "                  \"max_length\": 9.0,                             \n"
         "                  \"cutoff_radius\": 6.0                           \n"
         "                },                                                 \n"
         "                \"3\": {                                           \n"
         "                  \"max_length\": 8.0,                             \n"
         "                  \"cutoff_radius\": 6.0                           \n"
         "                },                                                 \n"
         "                ...                                                \n"
         "              }                                                  \n\n"

         "        orbit_specs: array (optional)                              \n"
         "          An array of clusters which are used to generate and      \n"
         "          include orbits of clusters whether or not they meet the  \n"
         "          `cutoff_radius` or `max_length` truncation criteria. See \n"
         "          the \"periodic_max_length\" docs for an example.\n\n";
}

/// Parse JSON to construct initial states for enumeration (as
/// std::vector<std::pair<std::string, ConfigEnumInput>>)
///
/// The result pairs a name, describing the Supercell or Configuration and
/// selected sites name, to a ConfigEnumInput.
///
/// This method enables several options for specifying initial states for
/// enumeration. There are two main categories of options: i) options specifying
/// supercells and configurations, and ii) options restricting which sites the
/// enumeration takes place on.
///
/// Specifying supercells and configurations: These options are all additive in
/// the sense that the result will be the combination of all that are found.
/// Selecting a Supercell as the initial state of enumeration is equivalent to
/// selecting the configuration in that supercell with all values of DoF set to
/// 0.
///
/// Options for specifying supercells and configurations are:
///     supercell_selection: string (optional)
///         Name of a selection of supercells to use as initial states for
///         enumeration
///     scelnames: array of string (optional)
///         Array of names of supercell to use as initial states for enumeration
///     config_selection: string (optional)
///         Name of a selection of configurations to use as initial states for
///         enumeration
///     confignames: array of string (optional)
///         Array of names of configuration to use as initial states for
///         enumeration
///     supercells: object (optional)
///         Specifies parameters for enumerating supercells. Options are:
///
///         min: int (optional, default=1)
///             Minimum volume supercells to enumerate
///         max: int (required)
///             Maximum volume supercells to enumerate
///         dirs: string (optional, default="abc")
///             Which lattice vectors of unit cell to enumerate over
///         unit_cell: 3x3 matrix of int or string (optional, default=identity
///         matrix)
///             The unit cell to tile into supercells. If string, use
///             transformation matrix for the supercell with given name.
///      fill/scelname, fill/supercell_selection, fill/supercells:
///          When enumerating local clusters using the "local_max_length"
///          method, prepend "background_" to specify the supercell(s) in which
///          to do the enumeration. Each supercell specified in this way is
///          filled with configurations specified using the other parameters.
///          Any configurations which do not exactly tile the specified
///          supercells are skipped.
///
/// Specifying all sites, particular sublattices, particular sites, or
/// particular clusters of sites to enumerate local DoF. The default behavior if
/// none of these options are given is selecting all sites for enumeration. The
/// values of DoF on all the sites that are not selected are frozen.
///
/// Options for selecting sites are:
///     sublats: array of integer (optional)
///         Indices of sublattices to allow enumeration on. May be used with
///         "sites".
///     sites: array of array of integer (optional)
///         Indices of sites to allow enumeration on, using [b, i, j, k]
///         notation (b=sublattice index, (i,j,k)=unit cell indices). May be
///         used with "sublats". Example:
///
///             "sites": [
///               [0, 0, 0, 0],
///               [0, 1, 0, 0],
///               [1, 0, 0, 0]
///             ]
///
///     cluster_specs: object (optional)
///         JSON object specifying orbits of clusters to generate. Each orbit
///         prototype is used to select sites to enumerate on each selected
///         supercell or configuration. If there are 4 supercells or
///         configurations selected, and there are 10 orbits generated, then
///         there will be 4*10=40 ConfigEnumInput generated. The "cluster_specs"
///         option cannot be used with the "sublats" or "sites" options. For
///         local clusters there may be even more generated, because of the
///         presence of distinct phenomenal clusters against the background
///         configuration.
///
///
/// Note:
/// - Parameter `PrimClex const *primclex` is used while transitioning Supercell
/// to no longer need
///   a `PrimClex const *`
/// - An error is inserted if no ConfigEnumInput are generated
void parse(
    InputParser<std::vector<std::pair<std::string, ConfigEnumInput>>> &parser,
    std::shared_ptr<Structure const> shared_prim, PrimClex const *primclex,
    DB::Database<Supercell> &supercell_db,
    DB::Database<Configuration> &configuration_db) {
  // check for "config_selection" and "confignames"
  DB::Selection<Configuration> config_selection;
  try {
    config_selection = DB::make_selection<Configuration>(
        configuration_db, parser.self, "confignames", "config_selection");
  } catch (std::exception &e) {
    std::stringstream msg;
    msg << "Error creating input states from configurations: " << e.what();
    parser.error.insert(msg.str());
  }

  // make_and_insert_canonical_supercell

  // option 1: create and use possibly non-canonical standalone supercell
  // (default) option 2: make_and_insert_canonical_supercell: bool (optional,
  // default = false)
  //           - make & insert canonical supercell, do not insert configuration
  // option 3: make_and_insert_canonical_configuration: true (optional, default
  // = false)
  //           - make & insert canonical supercell, and insert configuration

  // check for "config_list"
  typedef std::vector<std::pair<std::string, Configuration>> NamedConfiguration;
  NamedConfiguration config_list;

  if (parser.self.contains("config_list")) {
    try {
      auto it = parser.self.find("config_list");
      auto config_it = it->begin();
      auto config_end = it->end();
      for (; config_it != config_end; ++config_it) {
        auto config_ptr = config_it->make<Configuration>(shared_prim);
        std::string default_identifier = "unknown";
        if (primclex != nullptr) {
          config_ptr->supercell().set_primclex(primclex);
          default_identifier = config_ptr->name();
        }
        std::string identifier =
            config_it->get_if_else("identifier", default_identifier);
        config_list.emplace_back(identifier, std::move(*config_ptr));
      }
    } catch (std::exception &e) {
      std::stringstream msg;
      msg << "Error creating input states from config_list: " << e.what();
      parser.error.insert(msg.str());
    }
  }

  // check for "supercell_selection" and "scelnames"
  DB::Selection<Supercell> supercell_selection;
  try {
    supercell_selection = DB::make_selection<Supercell>(
        supercell_db, parser.self, "scelnames", "supercell_selection");
  } catch (std::exception &e) {
    std::stringstream msg;
    msg << "Error creating input states from supercells: " << e.what();
    parser.error.insert(msg.str());
  }
  // check for "supercells"
  auto scel_enum_props_subparser =
      parser.subparse_if<xtal::ScelEnumProps>("supercells", supercell_db);

  // check for "fill/"
  bool fill_supercells = false;
  DB::Selection<Supercell> fill_supercell_selection;
  std::shared_ptr<InputParser<xtal::ScelEnumProps>>
      fill_scel_enum_props_subparser;
  if (parser.self.contains("fill")) {
    fill_supercells = true;
    // check for "fill/supercell_selection" and "fill/scelnames"
    try {
      fill_supercell_selection =
          DB::make_selection<Supercell>(supercell_db, parser.self["fill"],
                                        "scelnames", "supercell_selection");
    } catch (std::exception &e) {
      std::stringstream msg;
      msg << "Error creating supercells to fill: " << e.what();
      parser.error.insert(msg.str());
    }

    // check for "fill/supercells"
    fill_scel_enum_props_subparser = parser.subparse_if<xtal::ScelEnumProps>(
        "fill/supercells", supercell_db);
  }

  // check for "sublats"
  std::vector<Index> sublats;
  parser.optional(sublats, "sublats");
  for (Index b : sublats) {
    if (b < 0 || b >= shared_prim->basis().size()) {
      std::stringstream msg;
      msg << "Error reading sublats: value out of range [0, "
          << shared_prim->basis().size() << ")";
      parser.error.insert(msg.str());
    }
  }

  // check for "sites"
  std::vector<UnitCellCoord> sites;
  parser.optional(sites, "sites");
  Index i = 0;
  for (UnitCellCoord site_uccoord : sites) {
    Index b = site_uccoord.sublattice();
    if (b < 0 || b >= shared_prim->basis().size()) {
      std::stringstream msg;
      msg << "Error reading sites[" << i
          << "]: sublattice index out of range [0, "
          << shared_prim->basis().size() << ")";
      parser.error.insert(msg.str());
    }
    ++i;
  }

  // Do not allow "cluster_specs" with "sublats" or "sites"
  //  (when would this be useful? which order to apply?)
  if ((sublats.size() || sites.size()) &&
      parser.self.contains("cluster_specs")) {
    std::stringstream msg;
    msg << "Error creating input states: "
        << "cannot include \"cluster_specs\" with \"sublats\" or \"sites\"";
    parser.error.insert(msg.str());
  }

  // check for "cluster_specs"
  auto cluster_specs_subparser = parser.subparse_if<ClusterSpecs>(
      "cluster_specs", shared_prim, shared_prim->factor_group());

  // at this point we have parsed everything except "cluster_specs",
  //   which is parsed later for each ConfigEnumInput
  if (!parser.valid()) {
    return;
  }

  // Use the supercells and configurations input to construct ConfigEnumInput,
  //   by default these have all sites selected
  std::vector<std::pair<std::string, ConfigEnumInput>> config_enum_input;
  for (const auto &config : config_selection.selected()) {
    config_enum_input.emplace_back(config.name(), config);
  }
  for (const auto &named_config : config_list) {
    config_enum_input.emplace_back(named_config.first, named_config.second);
  }
  for (const auto &scel : supercell_selection.selected()) {
    config_enum_input.emplace_back(scel.name(), scel);
  }
  if (scel_enum_props_subparser->value) {
    ScelEnumByProps enumerator{shared_prim, *scel_enum_props_subparser->value};
    for (auto const &supercell : enumerator) {
      auto supercell_it = supercell_db.insert(supercell).first;
      config_enum_input.emplace_back(supercell_it->name(), *supercell_it);
    }
  }

  // If "fill" is provided, the configurations of config_enum_input are used
  // to fill the provided supercells, skipping (with warning) any
  // configurations that do not exactly fill
  if (fill_supercells) {
    std::vector<std::pair<std::string, ConfigEnumInput>> tmp_config_enum_input;

    // lambda method fills supercell `scel` with all current
    // config_enum_input configurations
    auto _fill_supercell = [&](Supercell const &scel) {
      FillSupercell f{scel};
      for (auto const &name_config_pair : config_enum_input) {
        auto const &config = name_config_pair.second.configuration();
        SymOp const *op = f.find_symop(config.supercell().lattice());
        if (op == nullptr) {
          std::stringstream msg;
          msg << "Warning: input state " << name_config_pair.first
              << " cannot fill supercell " << scel.name();
          parser.warning.insert(msg.str());
          continue;
        }
        std::string name = name_config_pair.first;
        if (config.supercell().name() != scel.name()) {
          name = name + "-fill-" + scel.name();
        }
        tmp_config_enum_input.emplace_back(name, f(config));
      }
    };

    for (const auto &scel : fill_supercell_selection.selected()) {
      _fill_supercell(scel);
    }
    if (fill_scel_enum_props_subparser->value) {
      ScelEnumByProps enumerator{shared_prim,
                                 *fill_scel_enum_props_subparser->value};
      for (auto const &supercell : enumerator) {
        auto supercell_it = supercell_db.insert(supercell).first;
        _fill_supercell(*supercell_it);
      }
    }
    config_enum_input = std::move(tmp_config_enum_input);
  }
  // select sublattices and individual sites
  if (sublats.size() || sites.size()) {
    for (auto &input_name_value_pair : config_enum_input) {
      // select sublats & sites, and modify name with sublats & sites
      input_name_value_pair.second.clear_sites();

      std::stringstream ss;
      if (sublats.size()) {
        ss << ", sublats=" << jsonParser{sublats};
        input_name_value_pair.second.select_sublattices(sublats);
      }
      if (sites.size()) {
        ss << ", sites=" << jsonParser{sites};
        input_name_value_pair.second.select_sites(sites);
      }
      input_name_value_pair.first += ss.str();
    }
  }

  // select clusters
  if (cluster_specs_subparser->value != nullptr) {
    if (cluster_specs_subparser->value->periodicity_type() ==
        CLUSTER_PERIODICITY_TYPE::PRIM_PERIODIC) {
      // generate orbits from cluster_specs
      auto orbits =
          cluster_specs_subparser->value->make_periodic_orbits(CASM::log());

      // this will generate more ConfigEnumInput, with cluster sites selected
      std::vector<std::pair<std::string, ConfigEnumInput>>
          all_with_cluster_sites;
      for (auto &input_name_value_pair : config_enum_input) {
        std::vector<ConfigEnumInput> with_cluster_sites;
        select_cluster_sites(input_name_value_pair.second.configuration(),
                             orbits, std::back_inserter(with_cluster_sites));

        // modify name with cluster site indices
        for (auto const &value : with_cluster_sites) {
          std::stringstream ss;
          ss << ", cluster_sites=" << jsonParser{value.sites()};
          std::string name = input_name_value_pair.first + ss.str();
          all_with_cluster_sites.emplace_back(name, value);
        }
      }
      config_enum_input = std::move(all_with_cluster_sites);
    } else if (cluster_specs_subparser->value->periodicity_type() ==
               CLUSTER_PERIODICITY_TYPE::LOCAL) {
      std::cout << "WARNING: Local cluster enumeration is a work in progress"
                << std::endl;

      // TODO: generate inequivalent backgrounds
      bool include_all_inequivalent_reference_configs = true;
      // TODO: set phenomenal cluster site dof

      // generate orbits from cluster_specs
      auto orbits =
          cluster_specs_subparser->value->make_local_orbits(CASM::log());

      // this will generate more ConfigEnumInput, with cluster sites selected
      std::vector<std::pair<std::string, ConfigEnumInput>>
          all_with_cluster_sites;
      for (auto &input_name_value_pair : config_enum_input) {
        std::vector<ConfigEnumInput> with_cluster_sites;
        select_local_cluster_sites(input_name_value_pair.second.configuration(),
                                   orbits,
                                   include_all_inequivalent_reference_configs,
                                   std::back_inserter(with_cluster_sites));

        // modify name with cluster site indices
        for (auto const &value : with_cluster_sites) {
          std::stringstream ss;
          ss << ", cluster_sites=" << jsonParser{value.sites()};
          std::string name = input_name_value_pair.first + ss.str();
          all_with_cluster_sites.emplace_back(name, value);
        }
      }
      config_enum_input = std::move(all_with_cluster_sites);
    } else {
      std::stringstream msg;
      msg << "Error creating input states: unknown error using cluster_specs";
      parser.error.insert(msg.str());
    }
  }

  // Move all constructed ConfigEnumInput into the parser.value
  parser.value =
      notstd::make_unique<std::vector<std::pair<std::string, ConfigEnumInput>>>(
          std::move(config_enum_input));

  if (parser.value->size() == 0) {
    std::stringstream msg;
    msg << "Error creating input states: No supercells or configurations.";
    parser.error.insert(msg.str());
  }
}

}  // namespace CASM
