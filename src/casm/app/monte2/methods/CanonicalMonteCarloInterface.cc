#include "casm/app/monte2/methods/CanonicalMonteCarloInterface.hh"

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {

namespace Monte2 {

// // --- StateGenerator types ---
//
// enum class INITIAL_STATE_METHOD { DEFAULT, MIN_COST, BY_NAME, BY_CONFIG };
//
// struct InitialStateParams {
//   INITIAL_STATE_METHOD method;
//   std::function<double(GenericState)> calculator;
//   std::optional<std::string> configname;
//   std::optional<Configuration> configuration;
// };
//
// class StateGenerator {
//  public:
//   virtual GenericState next_state(
//       std::vector<ResultsSummary> const &results) const = 0;
//   virtual bool is_complete(
//       std::vector<ResultsSummary> const &results) const = 0;
// };
//
// class IncrementalStateGenerator : public StateGenerator {
//  public:
//   IncrementalCanonicalStateGenerator(
//       PrimClex const &primclex, InitialStateParams initial_state_params,
//       VectorValueMap begin_conditions, VectorValueMap increment_conditions,
//       VectorValueMap end_conditions, bool tol,
//       std::function<State(Configuration, VectorValueMap)> make_state,
//       bool dependent_runs);
//
//   // --- Begin required StateGenerator interface ---
//
//   State next_state(std::vector<ResultsSummary> const &results) const
//   override; bool is_complete(std::vector<ResultsSummary> const &results)
//   const override;
//
//   // --- End required StateGenerator interface ---
// };

/*
// --- Calculator types ---

class ClexCalculator {
 public:
  ClexCalculator(Clexulator clexulator, ECIContainer const &eci);

  // --- Begin required Calculator interface ---

  double value(State const &state) const;
  double delta_value(State const &state, Conversions const &convert,
                     OccEvent const &occ_event) const;

  // --- End required Calculator interface ---
};

class KerasCalculator {
 public:
  KerasCalculator(Clexulator clexulator,
                  std::vector<fs::path> asym_unit_models);

  // --- Begin required Calculator interface ---

  double value(State const &state) const;
  double delta_value(State const &state, Conversions const &convert,
                     OccEvent const &occ_event) const;

  // --- End required Calculator interface ---
};

class LammpsCalculator {
 public:
};
*/

namespace CanonicalMonteCarlo {

// // conditions:
// // "temperature_in_K"
// // "comp" (parametric composition) / "comp_n"
//
// // properties:
// // "extrinsic_formation_energy_in_eV"
//
// // default samplers:
// // - ... todo ...
//
// // eh
// bool check_conditions(VectorValueMap);
//
// // include PrimClex?
// ResultsSummary summarize(PrimClex const &primclex, Results const &results);
//
// std::vector<notstd::cloneable_ptr<Sampler>> default_samplers(
//     PrimClex const &primclex);
//
// struct Event {
//  public:
//   // Proposed occupation change
//   OccEvent occ_event;
//
//   // Change in formation energy due to occupation change
//   double delta_extrinsic_formation_energy_in_eV;
// };
//
// /// \brief Perform a Monte Carlo calculation at single point in thermodynamic
// /// conditions space.
// template <typename CalculatorType>
// Results single_run_canonical(State initial_state,
//                              CalculatorType formation_energy_calculator,
//                              CompletionCheck const &completion_check,
//                              Conversions convert,
//                              OccCandidateList occ_candidate_list, MTRand
//                              mtrand, SamplerMap sampler_map);
//
// /// Run Canonical Monte Carlo calculations
// template <typename CalculatorType>
// void run_canonical(CalculatorType formation_energy_calculator,
//                    StateGenerator const &state_generator,
//                    CompletionCheck const &completion_check, MTRand mtrand,
//                    SamplerMap sampler_map, ResultsIO const &results_io);

}  // namespace CanonicalMonteCarlo
}  // namespace Monte2

std::string CanonicalMonteCarloInterface::desc() const {
  // TODO
  std::string custom_options = "";

  std::string examples = "";

  return name() + ": \n\n" + custom_options + examples;
}

std::string CanonicalMonteCarloInterface::name() const {
  return "CanonicalMonteCarlo";
}

void CanonicalMonteCarloInterface::run(
    PrimClex &primclex, jsonParser const &json_options,
    jsonParser const &cli_options_as_json) const {
  Log &log = CASM::log();

  log.subsection().begin("CanonicalMonteCarlo");
  ParentInputParser parser =
      make_monte2_parent_parser(log, json_options, cli_options_as_json);
  std::runtime_error error_if_invalid{
      "Error reading CanonicalMonteCarlo JSON input"};

  // TODO

  log.custom("Checking input");

  log.custom("Checking for restart files");

  log.custom("Initializing Monte Carlo state");
}

// namespace Monte2 {
// namespace StandardCanonicalMonteCarlo {
//
// /// Run Canonical Monte Carlo calculations
// template <typename CalculatorType>
// void run_canonical(CalculatorType formation_energy_calculator,
//                    StateGenerator const &state_generator,
//                    CompletionCheck const &completion_check, MTRand mtrand,
//                    SamplerMap sampler_map, ResultsIO const &results_io) {
//   // enable restarts
//   std::vector<ResultsSummary> results_summary;
//   results_io.read_completed_runs_summary(results_summary);
//
//   // loop over thermodynamic conditions
//   while (!state_generator.is_complete(results_summary)) {
//     // get next initial state
//     State initial_state =
//         from_generic_state(state_generator.next_state(results_summary));
//     Supercell const &supercell = initial_state.configuration.supercell();
//     Conversions convert{supercell};
//     OccCandidateList occ_candidate_list{convert};
//
//     Results results = single_run_canonical(
//         initial_state, formation_energy_calculator, convert,
//         occ_candidate_list, mtrand, sampler_map, convergence_check);
//
//     // write results
//     results_summary.emplace_back(summarize(results));
//     results_io.write_completed_run(results, results_summary.size() - 1);
//   }
//
//   results_io.write_completed_runs_summary(results_summary);
// }
//
// }
// }

}  // namespace CASM
