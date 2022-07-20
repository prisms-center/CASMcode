#include "casm/app/sym/config_space_analysis.hh"

#include "casm/app/io/json_io.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/optional.hh"
#include "casm/clex/ConfigEnumByPermutation.hh"
#include "casm/clex/FillSupercell.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/io/json/Configuration_json_io.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/enumerator/DoFSpace.hh"
#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"

namespace CASM {
namespace DoFSpaceIO {

// included in dof_space_analysis.cc
void parse_dofs(ParentInputParser &parser, std::vector<DoFKey> &dofs,
                std::vector<DoFKey> const &all_dof_types);

}  // end namespace DoFSpaceIO

namespace ConfigSpaceIO {

jsonParser combine_config_space_json_options(
    jsonParser const &json_options, jsonParser const &cli_options_as_json) {
  std::map<std::string, std::string> cli_to_combined_keys{
      {"scelnames", "scelnames"},         // --scelnames
      {"confignames", "confignames"},     // --confignames
      {"selection", "config_selection"},  // --selection
      {"dofs", "dofs"}                    // --dofs
  };

  jsonParser json_combined{json_options};
  return combine_json_options(cli_to_combined_keys, cli_options_as_json,
                              json_combined);
}

Eigen::MatrixXd pseudoinverse(Eigen::MatrixXd const &M) {
  return M.transpose()
      .colPivHouseholderQr()
      .solve(Eigen::MatrixXd::Identity(M.cols(), M.cols()))
      .transpose();
}

}  // end namespace ConfigSpaceIO

/// Describe --config-space-analysis options
std::string config_space_analysis_desc() {
  std::string description =

      "Configuration space analysis (--config-space-analysis):            \n\n"

      "  The --config-space-analysis option (i) constructs a projector    \n"
      "  ont the DoF space spanned by all configurations symmetrically    \n"
      "  equivalent to a set of input configurations and (ii) finds the   \n"
      "  eigenvectors spanning that space. The eigenvectors are axes that \n"
      "  can be used as a basis for symmetry adapted order parameters.    \n"
      "  This method is faster than --dof-space-analysis, and often the   \n"
      "  resulting axes do lie along high-symmetry directions, but they   \n"
      "  may not lie within a single irreducible subspace, and they are   \n"
      "  not explicitly rotated to the highest symmetry directions.       \n\n";

  std::string custom_options =

      "  JSON options (--input or --settings):                            \n\n"

      "    dofs: string (optional, override with --dofs)                  \n"
      "      Names of degree of freedoms for which the analysis is run.   \n"
      "      The default includes all DoF types in the prim.              \n\n"

      "    tol: number (optional, default=1e-5)                           \n"
      "      Tolerance used for calculations.                             \n\n"

      "    exclude_homogeneous_modes: bool (optional, default=null)       \n"
      "      Exclude homogeneous modes if this is true, or include if     \n"
      "      this is false. If this is null (default), exclude homogeneous\n"
      "      modes for dof==\"disp\" only.                                \n\n"

      "    include_default_occ_modes: bool (optional, default=false)      \n"
      "      Include the dof component for the default occupation value on\n"
      "      each site with occupation DoF. The default is to exclude     \n"
      "      these modes because they are not independent. This parameter \n"
      "      is only checked dof==\"occ\".                                \n\n"

      "    output: string (optional, default=null)                        \n"
      "      If specified, output results to a file. Otherwise, print     \n"
      "      output.                                                      \n\n";

  return description + custom_options + parse_ConfigEnumInput_desc();
}

struct ConfigSpaceAnalysisResults {
  ConfigSpaceAnalysisResults(DoFSpace const &_standard_dof_space)
      : standard_dof_space(_standard_dof_space),
        projector(Eigen::MatrixXd::Zero(standard_dof_space.basis().cols(),
                                        standard_dof_space.basis().cols())) {}

  /// \brief Standard DoF space, may exclude default occupation modes
  ///     or homogeneous displacement modes, depending on method options
  DoFSpace standard_dof_space;

  /// \brief DoF values of all equivalent configurations in the
  ///     fully commensurate supercell, expressed in the basis of the
  ///     standard DoF space, with key == input configuration identifier
  std::map<std::string, std::vector<Eigen::VectorXd>> equivalent_dof_values;

  /// \brief All equivalent configurations in the fully commensurate
  ///     supercell, with key == input configuration identifier
  std::map<std::string, std::vector<Configuration>> equivalent_configurations;

  /// \brief Projection matrix
  Eigen::MatrixXd projector;

  /// \brief Non-zero eigenvalues of projector
  Eigen::VectorXd eigenvalues;

  /// \brief Symmetry-adapted config space, with basis formed by
  ///     eigenvectors of P
  std::unique_ptr<DoFSpace> symmetry_adapted_config_space;
};

/// \param dofs Types of DoF to run analysis for
/// \param
std::map<DoFKey, ConfigSpaceAnalysisResults> config_space_analysis(
    std::vector<DoFKey> dofs,
    std::map<std::string, Configuration> const &configurations,
    std::optional<bool> exclude_homogeneous_modes = std::nullopt,
    bool include_default_occ_modes = false, double tol = TOL) {
  std::map<DoFKey, ConfigSpaceAnalysisResults> results;

  if (configurations.size() == 0) {
    return results;
  }
  auto shared_prim = configurations.begin()->second.supercell().shared_prim();

  // --- Generate primitive configurations ---
  // prim config -> ID (might be duplicates from input configurations)
  std::map<Configuration, std::string> prim_configs;
  for (auto const &config : configurations) {
    prim_configs.emplace(config.second.primitive().in_canonical_supercell(),
                         config.first);
  }

  // --- Generate fully commensurate supercell ---
  std::set<xtal::Lattice> lattices;
  for (auto const &config : prim_configs) {
    lattices.insert(config.first.supercell().lattice());
  }
  auto const &fg = shared_prim->factor_group();
  xtal::Lattice super_lat = xtal::make_fully_commensurate_superduperlattice(
      lattices.begin(), lattices.end(), fg.begin(), fg.end());
  auto const &pg = shared_prim->point_group();
  super_lat = xtal::canonical::equivalent(super_lat, pg);
  auto shared_supercell =
      std::make_shared<Supercell const>(shared_prim, super_lat);

  // --- Generate symmetry adapted config spaces ---
  for (auto const &dof_key : dofs) {
    // --- Construct the standard DoF space ---
    DoFSpace dof_space(
        shared_prim, dof_key,
        shared_supercell->sym_info().transformation_matrix_to_super());

    if (!exclude_homogeneous_modes.has_value()) {
      if (dof_space.dof_key() == "disp") {
        exclude_homogeneous_modes = true;
      } else {
        exclude_homogeneous_modes = false;
      }
    }
    if (*exclude_homogeneous_modes) {
      dof_space = exclude_homogeneous_mode_space(dof_space);
    }

    if (dof_space.dof_key() == "occ" && !include_default_occ_modes) {
      dof_space = exclude_default_occ_modes(dof_space);
    }

    // --- Begin projector construction ---
    results.emplace(dof_key, dof_space);
    DoFSpace const &standard_dof_space = results.at(dof_key).standard_dof_space;
    auto &equivalent_dof_values = results.at(dof_key).equivalent_dof_values;
    auto &equivalent_configurations =
        results.at(dof_key).equivalent_configurations;
    Eigen::MatrixXd &P = results.at(dof_key).projector;

    // std::cout << "P:\n" << P << std::endl;
    for (auto const &prim_config : prim_configs) {
      Configuration prototype =
          fill_supercell(prim_config.first, shared_supercell);
      ConfigEnumByPermutation enumerator{prototype};
      std::vector<Eigen::VectorXd> equiv_x;
      std::vector<Configuration> equiv_config;
      for (auto const &config : enumerator) {
        equiv_config.push_back(config);
        Eigen::VectorXd x = get_normal_coordinate(config, standard_dof_space);

        // clean up x?
        for (int i = 0; i < x.size(); ++i) {
          if (almost_zero(x(i), tol)) {
            x(i) = 0.0;
          }
        }

        P += x * x.transpose();
        equiv_x.push_back(x);
        // std::cout << "P:\n" << P << std::endl;
      }
      equivalent_dof_values[prim_config.second] = equiv_x;
      equivalent_configurations[prim_config.second] = equiv_config;
    }
    // std::cout << "P:\n" << P << std::endl;

    // clean up P?
    for (int i = 0; i < P.rows(); ++i) {
      for (int j = 0; j < P.cols(); ++j) {
        if (almost_zero(P(i, j), tol)) {
          P(i, j) = 0.0;
        }
      }
    }

    // --- Eigendecomposition of P ---
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(P);
    Eigen::VectorXd D = solver.eigenvalues();
    Eigen::MatrixXd V = solver.eigenvectors();

    // std::cout << "eigenvalues:\n" << D << std::endl;
    // std::cout << "eigenvectors:\n" << V << std::endl;
    // std::cout << "check:\n" << V * D.asDiagonal() * V.transpose() <<
    // std::endl;

    // --- Identify non-zero eigenvalues and corresponding eigenvectors ---
    Eigen::VectorXd D_nonzero(D.size());
    Eigen::MatrixXd V_nonzero(P.rows(), P.cols());

    int i_nonzero = 0;
    for (int i = 0; i < D.size(); ++i) {
      if (!almost_zero(D(i), tol)) {
        D_nonzero(i_nonzero) = D(i);
        V_nonzero.col(i_nonzero) = V.col(i);
        ++i_nonzero;
      }
    }

    // std::cout << "non-zero eigenvalues:\n"
    //           << D_nonzero.head(i_nonzero) << std::endl;
    // std::cout << "non-zero eigenvectors:\n"
    //           << V_nonzero.leftCols(i_nonzero) << std::endl;
    // std::cout << "B * non-zero eigenvectors:\n"
    //           << standard_dof_space.basis() * V_nonzero.leftCols(i_nonzero)
    //           << std::endl;

    // --- Store results ---
    results.at(dof_key).eigenvalues = D_nonzero.head(i_nonzero);

    std::unique_ptr<DoFSpace> config_space;
    if (i_nonzero != 0) {
      results.at(dof_key).symmetry_adapted_config_space =
          std::make_unique<DoFSpace>(
              dof_space.shared_prim(), dof_space.dof_key(),
              dof_space.transformation_matrix_to_super(), dof_space.sites(),
              standard_dof_space.basis() * V_nonzero.leftCols(i_nonzero));
    }
  }

  return results;
}

/// Run config space analysis
///
void config_space_analysis(PrimClex &primclex, jsonParser const &json_options,
                           jsonParser const &cli_options_as_json) {
  using namespace DoFSpaceIO;
  using namespace ConfigSpaceIO;

  Log &log = CASM::log();
  std::shared_ptr<Structure const> shared_prim = primclex.shared_prim();

  log.subsection().begin<Log::debug>("config_space_analysis");
  log.indent() << "json_options:\n" << json_options << std::endl << std::endl;
  log.indent() << "cli_options_as_json:\n"
               << cli_options_as_json << std::endl
               << std::endl;
  log.end_section();

  // combine JSON options and CLI options
  jsonParser json_combined =
      combine_config_space_json_options(json_options, cli_options_as_json);

  // Read input data from JSON
  ParentInputParser parser{json_combined};
  std::runtime_error error_if_invalid{
      "Error reading `casm sym --config-space-analysis` input"};

  // 1) parse options

  // parse "dofs" (optional, default = all dof types)
  std::vector<std::string> dofs;
  parse_dofs(parser, dofs, all_dof_types(shared_prim->structure()));

  // parse "tol" (optional, default=1e-5)
  double tol = TOL;
  parser.optional(tol, "tol");

  // parse "output_dir" (optional, default = current_path)
  std::optional<fs::path> output;
  parser.optional(output, "output");

  // parse "exclude_homogeneous_modes" (optional, default = std::nullopt)
  std::optional<bool> exclude_homogeneous_modes;
  parser.optional(exclude_homogeneous_modes, "exclude_homogeneous_modes");

  // parse "include_default_occ_modes" (optional, default = false)
  bool include_default_occ_modes;
  parser.optional_else(include_default_occ_modes, "include_default_occ_modes",
                       false);

  // 2) parse input states
  typedef std::vector<std::pair<std::string, ConfigEnumInput>>
      NamedConfigEnumInput;
  auto input_parser_ptr = parser.parse_as<NamedConfigEnumInput>(
      shared_prim, &primclex, primclex.db<Supercell>(),
      primclex.db<Configuration>());
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  auto const &named_inputs = *input_parser_ptr->value;

  // 3) Run config space analysis
  std::map<std::string, Configuration> configurations;
  for (auto const &state : named_inputs) {
    configurations.emplace(state.first, state.second.configuration());
  }

  auto results =
      config_space_analysis(dofs, configurations, exclude_homogeneous_modes,
                            include_default_occ_modes, tol);

  // --- Report results ---
  jsonParser json;

  // put "config_list"
  json["config_list"].put_array();
  for (auto const &state : named_inputs) {
    jsonParser tjson;
    to_json(state.second.configuration(), tjson);
    tjson["identifier"] = state.first;
    json["config_list"].push_back(tjson);
  }

  // put results for each dof:
  for (auto const &val : results) {
    auto const &dof_key = val.first;
    auto const &result = val.second;

    // put "standard_dof_space"
    to_json(result.standard_dof_space, json[dof_key]["standard_dof_space"]);

    // put "equivalent_dof_values" && "all_dof_values"
    json[dof_key]["equivalent_dof_values"].put_obj();
    jsonParser &_all = json[dof_key]["all_dof_values"];
    _all.put_array();
    Index i_all = 0;
    for (auto const &pair : result.equivalent_dof_values) {
      auto const &identifier = pair.first;
      auto const &curr_dof_values = pair.second;
      jsonParser &_curr = json[dof_key]["equivalent_dof_values"][identifier];
      _curr.put_array();

      Index i = 0;
      for (auto const &dof_value : curr_dof_values) {
        _curr.push_back(jsonParser::null());
        to_json(dof_value, _curr[i], jsonParser::as_array());
        ++i;

        _all.push_back(jsonParser::null());
        to_json(dof_value, _all[i_all], jsonParser::as_array());
        ++i_all;
      }
    }

    // put "equivalent_configurations"
    json[dof_key]["equivalent_configurations"].put_obj();
    for (auto const &pair : result.equivalent_configurations) {
      auto const &identifier = pair.first;
      auto const &curr_configs = pair.second;
      jsonParser &_curr =
          json[dof_key]["equivalent_configurations"][identifier];
      _curr.put_array();

      Index i = 0;
      for (auto const &config : curr_configs) {
        _curr.push_back(jsonParser::object());
        to_json(config, _curr[i]);
        ++i;
      }
    }

    // put "projector"
    to_json(result.projector, json[dof_key]["projector"]);

    // put "eigenvalues"
    to_json(result.eigenvalues, json[dof_key]["eigenvalues"],
            jsonParser::as_array());

    // put "symmetry_adapted_config_space"
    if (result.symmetry_adapted_config_space == nullptr) {
      json[dof_key]["error"] = "Empty configuration space";
      json[dof_key]["symmetry_adapted_config_space"].put_null();
    } else {
      to_json(*result.symmetry_adapted_config_space,
              json[dof_key]["symmetry_adapted_config_space"]);
    }
  }

  // write results
  if (output.has_value()) {
    fs::ofstream sout(*output);
    sout << json << std::endl;
  } else {
    std::cout << json << std::endl;
  }
}

}  // namespace CASM
