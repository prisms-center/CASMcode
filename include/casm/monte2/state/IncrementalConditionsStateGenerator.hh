#ifndef CASM_monte2_IncrementalConditionsStateGenerator
#define CASM_monte2_IncrementalConditionsStateGenerator

#include <map>
#include <string>

#include "casm/monte2/definitions.hh"
#include "casm/monte2/state/State.hh"
#include "casm/monte2/state/StateGenerator.hh"

namespace CASM {
namespace Monte2 {

/// \brief Generates a series of states by constant conditions increments
///
/// The run information needed to check completion and generate subsequent
/// states is the vector of final states from the previous runs.
template <typename _ConfigType>
class IncrementalConditionsStateGenerator
    : public StateGenerator<_ConfigType, State<_ConfigType>> {
 public:
  using _ConfigType ConfigType;
  using State<_ConfigType> RunInfoType;

  /// \brief Constructor
  ///
  /// \param _initial_configuration The configuration for the initial state.
  /// \param _initial_conditions The conditions for the initial state.
  /// \param _conditions_increment The amount to change each conditions vector
  ///     between subsequent states.
  /// \param _n_states The total number of states to generate. Includes the
  ///     initial state.
  /// \param _dependent_runs If true, use the last configuration as the starting
  ///     point for the next state. If false, always use the configuration of
  ///     the initial state.
  IncrementalConditionsStateGenerator(
      State<ConfigType> const &_initial_state,
      VectorValueMap const &_conditions_increment, Index _n_states,
      bool _dependent_runs)
      : m_initial_state(_initial_state),
        m_conditions_increment(_conditions_increment),
        m_n_states(_n_states),
        m_dependent_runs(_dependent_runs) {}

  /// \brief Check if all requested states have been run
  bool is_complete(
      std::vector<State<ConfigType>> const &final_states) override {
    return final_states.size() >= m_n_states;
  }

  /// \brief Return the next state
  State<ConfigType> next_state(
      std::vector<State<ConfigType>> const &final_states) override {
    State<ConfigType> state((m_dependent_runs && final_states.size())
                                ? final_states.back()
                                : m_initial_state);
    if (final_states.size() == 0) {
      return state;
    } else {
      for (auto const &pair : m_conditions_increment) {
        state.conditions[pair.first] += pair.second;
      }
    }
    return state;
  }

 private:
  ConfigType m_initial_configuration;
  VectorValueMap m_initial_conditions;
  VectorValueMap m_conditions_increment;
  Index m_n_states;
  bool m_dependent_runs;
}

}  // namespace Monte2
}  // namespace CASM

#endif
