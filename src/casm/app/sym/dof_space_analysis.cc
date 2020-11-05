#include <map>
#include "casm/app/APICommand.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/io/json_io.hh"
#include "casm/app/sym/dof_space_analysis.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/enumerator/DoFSpace_impl.hh"
#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"
#include "casm/symmetry/SymRepTools.hh"
#include "casm/symmetry/io/json/SymRepTools.hh"

namespace dof_space_analysis_impl {

  using namespace CASM;

  jsonParser combine_dof_space_analysis_json_options(jsonParser const &json_options,
                                                     jsonParser const &cli_options_as_json) {

    std::map<std::string, std::string> cli_to_combined_keys {
      {"scelnames", "scelnames"},         // --scelnames
      {"confignames", "confignames"},     // --confignames
      {"selection", "config_selection"},  // --selection
      {"dofs", "dofs"},                   // --dofs
      {"calc_wedge", "calc_wedge"}        // --calc-wedge
    };

    jsonParser json_combined {json_options};
    return combine_json_options(cli_to_combined_keys,
                                cli_options_as_json,
                                json_combined);
  }

  template<typename PermuteIteratorIt>
  void write_config_symmetry_files(ConfigEnumInput const &config,
                                   PermuteIteratorIt group_begin,
                                   PermuteIteratorIt group_end,
                                   fs::path sym_dir) {

    Lattice config_lattice = config.configuration().ideal_lattice();

    SymGroup config_factor_group = make_sym_group(group_begin, group_end, config_lattice);

    // Write lattice point group
    {
      SymGroup config_lattice_pg(SymGroup::lattice_point_group(config_lattice));
      jsonParser json;
      write_symgroup(config_lattice_pg, json);
      fs::ofstream outfile {sym_dir / "lattice_point_group.json"};
      json.print(outfile);
    }

    // Write factor group
    {
      jsonParser json;
      write_symgroup(config_factor_group, json);
      fs::ofstream outfile {sym_dir / "factor_group.json"};
      json.print(outfile);
    }

    // Write crystal point group
    {
      SymGroup config_point_group = make_point_group(group_begin, group_end, config_lattice);
      jsonParser json;
      write_symgroup(config_point_group, json);
      fs::ofstream outfile {sym_dir / "crystal_point_group.json"};
      json.print(outfile);
    }
  }

  /// For now, only support "confignames" and "config_selection" for reading ConfigEnumInput
  ///
  /// Later, could support other ConfigEnumInput JSON options (supercells, selecting sites, clusters,
  /// etc.), but need to determine how to write / print the data
  void require_database_configurations(ParentInputParser &parser) {
    std::set<std::string> do_not_allow {
      "scelnames",
      "supercell_selection",
      "supercells",
      "sublats",
      "sites",
      "cluster_specs"};

    for(auto key : do_not_allow) {
      if(parser.self.contains(key)) {
        std::stringstream msg;
        msg << "Error in dof_space_analysis: \"" << key << "\" is not an allowed option.";
        parser.error.insert(msg.str());
      }
    }
  }

  /// Parser "dofs" value for dof space analysis
  ///
  /// dofs: array of string (optional, default=all_dof_types)
  ///     Entries must exist in "all_dof_types" else an error is inserted.
  ///
  void parse_dofs(ParentInputParser &parser,
                  std::vector<DoFKey> &dofs,
                  std::vector<DoFKey> const &all_dof_types) {
    parser.optional_else(dofs, "dofs", all_dof_types);
    for(DoFKey const &dof : dofs) {
      if(std::find(all_dof_types.begin(), all_dof_types.end(), dof) == all_dof_types.end()) {
        std::stringstream msg;
        msg << "Error parsing \"dofs\": \"" << dof << "\" is invalid.";
        parser.error.insert(msg.str());
      }
    }
  }

}

namespace CASM {

  /// Describe DoF space analysis
  std::string dof_space_analysis_desc() {
    std::string description =

      "Analyze DoF spaces (--dof-space-analysis): \n\n"

      "  The --dof-space-analysis option writes a symmetry analysis of one or more     \n"
      "  DoF spaces.                                                                   \n\n"

      "  Usage:                                                                      \n"
      "  - Specify one or more DoF spaces and a symmetry group via selection of      \n"
      "    configurations and DoF type.                                              \n"
      "    - The supercell, selected sites, DoF type, and subspace determine the     \n"
      "      DoF space.                                                              \n"
      "    - The configuration factor group (symmetry operations that keep the       \n"
      "      supercell lattice and the configuration DoF values invariant) is used   \n"
      "      for symmetry analysis.                                                  \n"
      "  - The analysis finds irreducible representations of the DoF vector space and\n"
      "    optionally calculates the symmetrically unique wedges in that space.      \n\n"

      "  JSON options (--input or --settings):                                       \n\n"

      "    dofs: array of string (optional, override with --dofs)                       \n"
      "      Name of degree of freedoms for which the DoF space analysis is performed   \n"
      "      for each input configuration. The default includes all DoF types in the    \n"
      "      prim.                                                                      \n\n"

      "    confignames: Array of strings (optional, override with --confignames)        \n"
      "      Names of configurations to be used as initial states for enumeration. All  \n"
      "      specified sublattices or sites will be enumerated on and all other DoFs    \n"
      "      will maintain the values of the initial state.                             \n"
      "      Ex: \"confignames\" : [\"SCEL1_1_1_1_0_0_0/1\",\"SCEL2_2_1_1_0_0_0/3\"]    \n\n"

      "    config_selection: string (optional, override with --selection)               \n"
      "      Name of a selection of configurations to perform analysis on.              \n\n"

      "    calc_wedge: bool (optional, default=false, override with --calc-wedge)       \n"
      "      Perform calculation of irreducible wedge (may significantly slow down      \n"
      "      analysis).                                                                 \n\n";

    return description;
  }

  /// Perform DoF space analysis
  ///
  /// Usage:
  /// - Specify one or more DoF spaces and a symmetry group via selection of configurations,
  ///   sites within the configurations, DoF type, and DoF subspace.
  ///   - The supercell, selected sites, DoF type, and subspace determine the DoF space
  ///   - The subgroup of the configuration factor group that leaves selected sites invariant
  ///     is used for symmetry analysis
  /// - The analysis finds irreducible representations of the DoF vector space and
  ///   optionally calculates the symmetrically unique wedges in that space
  /// - See `dof_space_analysis_desc()` for JSON and CLI input options
  void dof_space_analysis(PrimClex &primclex,
                          jsonParser const &json_options,
                          jsonParser const &cli_options_as_json) {

    using namespace dof_space_analysis_impl;

    Log &log = CASM::log();
    DirectoryStructure const &dir = primclex.dir();

    log.subsection().begin<Log::debug>("dof_space_analysis");
    log.indent() << "json_options:\n" << json_options << std::endl << std::endl;
    log.indent() << "cli_options_as_json:\n" << cli_options_as_json << std::endl << std::endl;
    log.end_section();

    // combine JSON options and CLI options
    jsonParser json_combined = combine_dof_space_analysis_json_options(
                                 json_options,
                                 cli_options_as_json);

    log.indent() << "Input:\n" << json_combined << std::endl << std::endl;

    // Read input data from JSON
    ParentInputParser parser {json_combined};
    std::runtime_error error_if_invalid {"Error reading `casm sym --dof-space-analysis` input"};

    // For now, only allow input that is configurations already existing in the database
    require_database_configurations(parser);

    typedef std::vector<std::pair<std::string, ConfigEnumInput>> NamedConfigEnumInput;
    auto input_parser_ptr = parser.parse_as<NamedConfigEnumInput>(
                              primclex.shared_prim(),
                              &primclex,
                              primclex.db<Supercell>(),
                              primclex.db<Configuration>());

    // parse "dofs" (optional, default = all dof types)
    std::vector<DoFKey> dofs;
    parse_dofs(parser, dofs, all_dof_types(primclex.prim().structure()));

    // parse "calc_wedge" (optional, default = false)
    bool calc_wedge;
    parser.optional_else(calc_wedge, "calc_wedge", false);

    report_and_throw_if_invalid(parser, log, error_if_invalid);

    // For each enumeration envrionment, for each DoF type specified, perform analysis and write files.
    auto const &named_inputs = *input_parser_ptr->value;
    for(auto const &named_input : named_inputs) {

      std::string name = named_input.first;
      log.begin(name);
      log.increase_indent();

      ConfigEnumInput const &config_enum_input = named_input.second;
      Configuration const &configuration = config_enum_input.configuration();
      std::vector<PermuteIterator> group = make_invariant_subgroup(configuration);

      // These should not occur for now. They should be prevented by require_database_configurations.
      // TODO: add support for other dof space analyses
      if(configuration.id() == "none") {
        throw std::runtime_error("Error in dof_space_analysis: configuration does not exist in database.");
      }
      if(name != configuration.name()) {
        throw std::runtime_error("Error in dof_space_analysis: name error.");
      }
      if(config_enum_input.sites().size() != configuration.size()) {
        throw std::runtime_error("Error in dof_space_analysis: incomplete site selection error.");
      }

      fs::path sym_dir = dir.symmetry_dir(configuration.name());
      fs::create_directories(sym_dir);

      // for configuration, write lattice_point_group.json, factor_group.json, crystal_point_group.json
      write_config_symmetry_files(config_enum_input, group.begin(), group.end(), sym_dir);

      // write "dof_analysis_<dof>.json" for each specified DoF type
      for(DoFKey const &dof : dofs) {
        log << "Working on: " << name << " " << dof << std::endl;
        DoFSpace dof_space {config_enum_input, dof};
        auto report = vector_space_sym_report(dof_space, group.begin(), group.end(), calc_wedge);

        jsonParser json;
        to_json(report, json);
        to_json(dof_space, json, name);

        std::string filename = "dof_analysis_" + dof + ".json";
        json.write(sym_dir / filename);
        log << "Writing: " << (sym_dir / filename) << std::endl << std::endl;

      }

      log.decrease_indent();
      log << std::endl;
    }
  }
}
