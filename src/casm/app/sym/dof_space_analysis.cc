#include "casm/app/sym/dof_space_analysis.hh"

#include "casm/app/io/json_io.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/optional.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/enumerator/DoFSpace.hh"
#include "casm/enumerator/io/dof_space_analysis.hh"
#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"

namespace CASM {
namespace DoFSpaceIO {

jsonParser combine_dof_space_analysis_json_options(
    jsonParser const &json_options, jsonParser const &cli_options_as_json) {
  std::map<std::string, std::string> cli_to_combined_keys{
      {"scelnames", "scelnames"},         // --scelnames
      {"confignames", "confignames"},     // --confignames
      {"selection", "config_selection"},  // --selection
      {"dofs", "dofs"},                   // --dofs
      {"calc_wedge", "calc_wedge"}        // --calc-wedge
  };

  jsonParser json_combined{json_options};
  return combine_json_options(cli_to_combined_keys, cli_options_as_json,
                              json_combined);
}

/// For output_type==symmetry_directory, the input is restricted to
/// configurations currently existing in the database
void require_database_configurations(ParentInputParser &parser) {
  std::set<std::string> do_not_allow{"scelnames",  "supercell_selection",
                                     "supercells", "sublats",
                                     "sites",      "cluster_specs"};

  for (auto key : do_not_allow) {
    if (parser.self.contains(key)) {
      std::stringstream msg;
      msg << "Error in dof_space_analysis: \"" << key
          << "\" is not an allowed option when 'output_type' == "
             "'symmetry_directory' (only configurations currently existing in "
             "the database are allowed).";
      parser.error.insert(msg.str());
    }
  }
}

/// Parser "dofs" value for dof space analysis
///
/// dofs: array of string (optional, default=all_dof_types)
///     Entries must exist in "all_dof_types" else an error is inserted.
///
void parse_dofs(ParentInputParser &parser, std::vector<DoFKey> &dofs,
                std::vector<DoFKey> const &all_dof_types) {
  parser.optional_else(dofs, "dofs", all_dof_types);
  for (DoFKey const &dof : dofs) {
    if (std::find(all_dof_types.begin(), all_dof_types.end(), dof) ==
        all_dof_types.end()) {
      std::stringstream msg;
      msg << "Error parsing \"dofs\": \"" << dof << "\" is invalid.";
      parser.error.insert(msg.str());
    }
  }
}

}  // end namespace DoFSpaceIO

/// Describe DoF space analysis
std::string dof_space_analysis_desc() {
  std::string description =

      "Analyze DoF spaces (--dof-space-analysis): \n\n"

      "  The --dof-space-analysis option writes a symmetry analysis of one or "
      "more     \n"
      "  DoF spaces.                                                           "
      "        \n\n"

      "  Usage:                                                                "
      "      \n"
      "  - Specify one or more DoF spaces and a symmetry group via selection "
      "of      \n"
      "    configurations and DoF type.                                        "
      "      \n"
      "    - The supercell, selected sites, DoF type, and subspace determine "
      "the     \n"
      "      DoF space.                                                        "
      "      \n"
      "    - The configuration factor group (symmetry operations that keep the "
      "      \n"
      "      supercell lattice and the configuration DoF values invariant) is "
      "used   \n"
      "      for symmetry analysis.                                            "
      "      \n"
      "  - The analysis finds irreducible representations of the DoF vector "
      "space and\n"
      "    optionally calculates the symmetrically unique wedges in that "
      "space.      \n\n";

  std::string custom_options =

      "  JSON options (--input or --settings):                                 "
      "      \n\n"

      "    dofs: array of string (optional, override with --dofs)              "
      "         \n"
      "      Name of degree of freedoms for which the DoF space analysis is "
      "performed   \n"
      "      for each input configuration. The default includes all DoF types "
      "in the    \n"
      "      prim.                                                             "
      "         \n\n"

      "    calc_wedge: bool (optional, default=false, override with "
      "--calc-wedge)       \n"
      "      Perform calculation of irreducible wedge (may significantly slow "
      "down      \n"
      "      analysis).                                                        "
      "         \n\n"

      "    exclude_homogeneous_modes: bool (optional, default=null)       \n"
      "      Exclude homogeneous modes if this is true, or include if     \n"
      "      this is false. If this is null (default), only exclude       \n"
      "      homogeneous modes for dof==\"disp\" and the input space      \n"
      "      includes all supercell sites.                                \n\n"

      "    include_default_occ_modes: bool (optional, default=false)      \n"
      "      Include the dof component for the default occupation value on\n"
      "      each site with occupation DoF. The default is to exclude     \n"
      "      these modes because they are not independent. This parameter \n"
      "      is only checked if DoF is \"occ\" and \"axes\" are not       \n"
      "      included explicitly.  \n\n"

      "    axes: matrix or JSON object (optional)                           \n"
      "      Coordinate axes of the DoF grid. This parameter is only checked\n"
      "      when there is a single input state and single \"dof\". The     \n"
      "      default value is the identity matrix of DoF space dimension.   \n"
      "      Each element in an axis vector correponds to an individual DoF.\n"
      "      Each vector corresponds to a collective mode. If not included, \n"
      "      the \full space is included and a glossary describing which DoF\n"
      "      is specified by which vector element is generated. The 'axes'  \n"
      "      may be rank deficient, specifying a subspace of the full DoF   \n"
      "      space specified by the \"dof\" value and initial state.      \n\n"

      "      Example if matrix (row vector matix): \n"
      "        \"axes\" : [ \n"
      "          [1, 1, 1, 1, 1, 1], \n"
      "          [1,-1, 0,-1, 1, 0], \n"
      "          [1,-1, 0, 1,-1, 0]  \n"
      "        ] \n\n"

      "      Example if JSON object (named axis vectors): \n"
      "          \"axes\": { \n"
      "            \"q1\": [1, 1, 1, 1, 1, 1], \n"
      "            \"q2\": [1,-1, 0,-1, 1, 0], \n"
      "            \"q3\": [1,-1, 0, 1,-1, 0]  \n"
      "          } \n\n"

      "      Note: \n"
      "      - If some \"qi\" in the range [1, DoF space dimension] are     \n"
      "        missing, then the subspace is generated using the axes       \n"
      "        that are provided.                                           "
      "\n\n"

      "    write_symmetry: bool (optional, default=true)                       "
      "         \n"
      "      If true (default), write the lattice point group, factor group "
      "(operations \n"
      "      leave the configuration and selected sites invariant), and "
      "crystal point   \n"
      "      group (factor group operations excluding translations) of the "
      "initial state\n"
      "      for analysis.                                                     "
      "         \n\n"

      "    write_structure: bool (optional, default=true)                      "
      "         \n"
      "      If true (default), write a \"structure.json\" file containing the "
      "structure\n"
      "      generated from the configuration DoF.                             "
      "         \n\n"

      "    output_type: string (optional, default=\"symmetry_directory\")      "
      "         \n"
      "      Selects how output files are written. Options are:                "
      "         \n\n"

      "      \"symmetry_directory\": (default)                            \n"
      "        If selected, only accepts \"confignames\" and              \n"
      "        \"config_selection\" to specify the initial state for      \n"
      "        analysis and the results are stored in dedicated folders in\n"
      "        the CASM project symmetry directory:                       \n"
      "            \"<project_path>/symmetry/analysis/<configname>\".     \n\n"

      "      \"sequential\":                                              \n"
      "        If selected, accept any input for specifying the initial   \n"
      "        state for analysis, including \"scelnames\",               \n"
      "        \"supercell_selection\", \"supercells\"                    \n"
      "        \"sublats\", \"sites\", and \"cluster_specs\", and the     \n"
      "        results are stored in indexed folders                      \n"
      "        \"<output_dir>/dof_space_analysis/state.<index>\".         \n\n"

      "      \"combined_json_file\":                                      \n"
      "        If selected, accept any input for specifying the initial   \n"
      "        state for analysis, including \"scelnames\",               \n"
      "        \"supercell_selection\", \"supercells\", \"sublats\",      \n"
      "        \"sites\", and \"cluster_specs\", and the results are      \n"
      "        stored in a single JSON file \"<output_dir>/dof_space.json\"\n"
      "        instead of being written separately.                       \n\n"

      "      \"combined_json\":                                           \n"
      "        If selected, accept any input for specifying the initial   \n"
      "        state for analysis, including \"scelnames\",               \n"
      "        \"supercell_selection\", \"supercells\", \"sublats\",      \n"
      "        \"sites\", and \"cluster_specs\", and the results are      \n"
      "        output as a single JSON array.                             \n\n"

      "    output_dir: string (optional, default=current path)            \n"
      "      Selects where output files are written.                      \n\n";

  return description + custom_options + parse_ConfigEnumInput_desc();
}

/// Perform DoF space analysis
///
/// Usage:
/// - Specify one or more DoF spaces and a symmetry group via selection of
/// configurations,
///   sites within the configurations, DoF type, and DoF subspace.
///   - The supercell, selected sites, DoF type, and subspace determine the DoF
///   space
///   - The subgroup of the configuration factor group that leaves selected
///   sites invariant
///     is used for symmetry analysis
/// - The analysis finds irreducible representations of the DoF vector space and
///   optionally calculates the symmetrically unique wedges in that space
/// - See `dof_space_analysis_desc()` for JSON and CLI input options
void dof_space_analysis(PrimClex &primclex, jsonParser const &json_options,
                        jsonParser const &cli_options_as_json) {
  using namespace DoFSpaceIO;

  Log &log = CASM::log();

  log.subsection().begin<Log::debug>("dof_space_analysis");
  log.indent() << "json_options:\n" << json_options << std::endl << std::endl;
  log.indent() << "cli_options_as_json:\n"
               << cli_options_as_json << std::endl
               << std::endl;
  log.end_section();

  // combine JSON options and CLI options
  jsonParser json_combined = combine_dof_space_analysis_json_options(
      json_options, cli_options_as_json);

  // Read input data from JSON
  ParentInputParser parser{json_combined};
  std::runtime_error error_if_invalid{
      "Error reading `casm sym --dof-space-analysis` input"};

  // parse "output_type"= "symmetry_directory" (default), "sequential", or
  // "combined_json_file", or "combined_json"
  std::string output_type;
  parser.optional_else(output_type, "output_type",
                       std::string("symmetry_directory"));

  if (output_type != "combined_json") {
    log.indent() << "Input:\n" << json_combined << std::endl << std::endl;
  }

  // 1) parse options

  DoFSpaceAnalysisOptions options;

  // parse "dofs" (optional, default = all dof types)
  parse_dofs(parser, options.dofs, all_dof_types(primclex.prim().structure()));

  // parse "calc_wedge" (optional, default = false)
  parser.optional_else(options.calc_wedge, "calc_wedge", false);

  // parse "write_symmetry" (optional, default = true)
  parser.optional_else(options.write_symmetry, "write_symmetry", true);

  // parse "write_structure" (optional, default = true)
  parser.optional_else(options.write_structure, "write_structure", true);

  // parse "output_dir" (optional, default = current_path)
  fs::path output_dir;
  parser.optional_else(output_dir, "output_dir", fs::current_path());

  // parse "exclude_homogeneous_modes" (optional, default = std::nullopt)
  parser.optional(options.exclude_homogeneous_modes,
                  "exclude_homogeneous_modes");

  // parse "include_default_occ_modes" (optional, default = false)
  parser.optional_else(options.include_default_occ_modes,
                       "include_default_occ_modes", false);

  // 2) parse input states

  typedef std::vector<std::pair<std::string, ConfigEnumInput>>
      NamedConfigEnumInput;
  auto input_parser_ptr = parser.parse_as<NamedConfigEnumInput>(
      primclex.shared_prim(), &primclex, primclex.db<Supercell>(),
      primclex.db<Configuration>());
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  auto const &named_inputs = *input_parser_ptr->value;

  // 3) parse "axes"
  if (named_inputs.size() == 1 && options.dofs.size() == 1 &&
      parser.self.contains("axes")) {
    ConfigEnumInput const &config_input = named_inputs[0].second;
    Supercell const &supercell = config_input.configuration().supercell();
    Index dof_space_dimension = get_dof_space_dimension(
        options.dofs[0], supercell.prim(),
        supercell.sym_info().transformation_matrix_to_super(),
        config_input.sites());

    jsonParser tmp = json_combined;
    tmp["min"] = 0.0;
    tmp["max"] = 1.0;
    tmp["num"] = 1;
    InputParser<AxesCounterParams> axes_parser(tmp, dof_space_dimension);
    report_and_throw_if_invalid(axes_parser, log, error_if_invalid);
    options.basis = axes_parser.value->axes;
  }

  // 4) Construct output method implementation and run dof space analysis

  if (output_type == "combined_json_file") {
    // write all output to one large JSON file:
    // <output_dir>/dof_space.json
    CombinedJsonOutput output{output_dir};
    dof_space_analysis(named_inputs, options, output);
  } else if (output_type == "combined_json") {
    // write all output to log as one JSON array
    CombinedJsonOutput output;
    dof_space_analysis(named_inputs, options, output);
    log << output.combined_json() << std::endl;
  } else if (output_type == "sequential") {
    // write output sequentially indexed directories:
    // <output_dir>/dof_space_analysis/state.<index>/
    SequentialDirectoryOutput output{output_dir};
    dof_space_analysis(named_inputs, options, output);
  } else if (output_type == "symmetry_directory") {
    require_database_configurations(parser);
    report_and_throw_if_invalid(parser, log, error_if_invalid);

    // write output to CASM project symmetry directory:
    // <project_path>/symmetry/analysis/<configname>/
    SymmetryDirectoryOutput output{primclex.dir()};
    dof_space_analysis(named_inputs, options, output);
  } else {
    std::string msg =
        "Error: 'output_type' must be one of \"symmetry_directory\" (default), "
        "\"sequential\", or \"combined_json\"";
    parser.insert_error("output_type", msg);
    report_and_throw_if_invalid(parser, log, error_if_invalid);
  }
}
}  // namespace CASM
