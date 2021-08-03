#ifndef CASM_monte2_methods_Canonical
#define CASM_monte2_methods_Canonical

namespace CASM {
namespace Monte2 {
namespace Canonical {

// conditions:
// "temperature_in_K"
// "comp" (parametric composition) / "comp_n"

// properties:
// "extrinsic_formation_energy_in_eV"

// default samplers:
// - ... todo ...

/// Check that input conditions have required names and sizes
///
/// Input:
/// - "temperature": required, shape=(1,)
///   - temperature in K
/// - "comp": optional
///   - Parametric composition. Must be consistent with composition axes.
///     One of "comp" or "comp_n" is required.
/// - "comp_n": optional
///   - Composition per unit cell. Must be consistent with composition axes.
///     One of "comp" or "comp_n" is required.
///
/// Standard form after validation:
/// - "temperature": shape=(1,)
///   - temperature in K
/// - "comp_n":
///   - shape determined from composition_axes
///
/// \params shared_prim Prim for this calculation
/// \params composition_axes Composition axes for this system
/// \params conditions Name-value pairs, see requirements
///
/// \returns A (shared) input validator, with any error or warning messages,
/// and, if valid, the conditions put into standard form
///
std::shared_ptr<InputValidator<VectorValueMap>> validate_conditions(
    CompositionConverter const &composition_axes,
    VectorValueMap const &conditions);

/// \brief Perform a Monte Carlo calculation at single point in thermodynamic
/// conditions space.
template <typename CalculatorType>
Results single_run_canonical(State initial_state,
                             CalculatorType formation_energy_calculator,
                             CompletionCheck completion_check,
                             Conversions convert,
                             OccCandidateList occ_candidate_list,
                             MTRand &mtrand,
                             StateSamplingFunctionMap sampling_functions);

/// Run Canonical Monte Carlo calculations
template <typename CalculatorType>
void run_canonical(CalculatorType formation_energy_calculator,
                   StateGenerator state_generator,
                   CompletionCheck completion_check, MTRand &mtrand,
                   StateSamplingFunctionMap sampling_functions,
                   ResultsIO results_io) {
  std::vector<Results> results_summary;

  // Checks if restarting
  results_io.read_results_summary(results_summary);

  // loop over thermodynamic conditions
  while (!state_generator.is_complete(results_summary)) {
    // get next initial state
    State initial_state = state_generator.next_state(results_summary);

    // Prepare:
    // - the index converter
    // - the list of occupants allowed to swap
    Supercell const &supercell = initial_state.configuration.supercell();
    Conversions convert{supercell};
    OccCandidateList occ_candidate_list{convert};

    // Run Monte Carlo at a single condition
    Results results = single_run_canonical(
        initial_state, formation_energy_calculator, convert, occ_candidate_list,
        mtrand, sampling_functions, convergence_check);

    // Write results for this condition
    results_io.write_results(results, results_summary.size() - 1);

    // Summarize (erase detailed sampled data) and write results summary
    results.sampled_data.reset();
    results_summary.emplace_back(results);
    results_io.write_results_summary(results_summary);
  }
}

}  // namespace Canonical
}  // namespace Monte2
}  // namespace CASM

#endif
