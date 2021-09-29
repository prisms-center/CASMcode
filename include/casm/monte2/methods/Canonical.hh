#ifndef CASM_monte2_methods_Canonical
#define CASM_monte2_methods_Canonical

namespace CASM {
namespace Monte2 {
namespace Canonical {

// /// Check that input conditions have required names and sizes
// ///
// /// Input:
// /// - "temperature": required, shape=(1,)
// ///   - temperature in K
// /// - "comp": optional
// ///   - Parametric composition. Must be consistent with composition axes.
// ///     One of "comp" or "comp_n" is required.
// /// - "comp_n": optional
// ///   - Composition per unit cell. Must be consistent with composition axes.
// ///     One of "comp" or "comp_n" is required.
// ///
// /// Standard form after validation:
// /// - "temperature": shape=(1,)
// ///   - temperature in K
// /// - "comp_n":
// ///   - shape determined from composition_axes
// ///
// /// \params shared_prim Prim for this calculation
// /// \params composition_axes Composition axes for this system
// /// \params conditions Name-value pairs, see requirements
// ///
// /// \returns A (shared) input validator, with any error or warning messages,
// /// and, if valid, the conditions put into standard form
// ///
// std::shared_ptr<InputValidator<VectorValueMap>> validate_conditions(
//     CompositionConverter const &composition_axes,
//     VectorValueMap const &conditions);

// /// Run Canonical Monte Carlo calculations
// template <typename CalculatorType>
// void run_canonical(CalculatorType formation_energy_calculator,
//                    StateGenerator state_generator,
//                    CompletionCheck completion_check, MTRand &mtrand,
//                    StateSamplingFunctionMap sampling_functions,
//                    ResultsIO results_io);

}  // namespace Canonical
}  // namespace Monte2
}  // namespace CASM

#endif
