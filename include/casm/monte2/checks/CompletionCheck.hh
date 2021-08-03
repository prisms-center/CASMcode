#ifndef CASM_monte2_CompletionCheck
#define CASM_monte2_CompletionCheck

#include <optional>

#include "casm/monte2/Definitions.hh"
#include "casm/monte2/Sampling.hh"
#include "casm/monte2/checks/ConvergenceCheck.hh"
#include "casm/monte2/checks/CutoffCheck.hh"
#include "casm/monte2/checks/EquilibrationCheck.hh"

namespace CASM {
namespace Monte2 {

// --- Completion checking (cutoff & convergence) ---

/// \brief Parameters that determine if a calculation is complete
struct CompletionCheckParams {
  /// \brief Completion check parameters that don't depend on the sampled values
  CutoffCheckParams cutoff_params;

  /// \brief Sampler components that must be checked for convergence
  std::map<SamplerComponent, ConvergenceCheckParams> convergence_check_params;

  /// \brief Minimum number of samples before checking for completion
  CountType min_sample = 0;

  /// \brief How often to check for completion
  ///
  /// Check for completion performed if:
  /// - n_samples % check_frequency == 0
  /// - and n_samples >= min_sample
  CountType check_frequency = 0;
};

/// \brief Stores completion check results
struct CompletionCheckResults {
  /// True if calculation is complete, either due to convergence or cutoff
  bool is_complete = false;

  EquilibrationCheckResults equilibration_check_results;

  ConvergenceCheckResults convergence_check_results;
};

/// \brief Checks if a cutoff or convergence criteria are met
class CompletionCheck {
 public:
  CompletionCheck(CompletionCheckParams params);

  // --- Begin required CompletionCheck interface ---

  bool is_complete(SampledData const &sampled_data);

  CompletionCheckResults const &results() const { return m_results; }

  // --- End required CompletionCheck interface ---

 private:
  void _check(SampledData const &sampled_data);

  CompletionCheckParams m_params;

  CompletionCheckResults m_results;
};

// --- Inline definitions ---

inline bool CompletionCheck::is_complete(SampledData const &sampled_data) {
  CountType n_samples = get_n_samples(sampled_data);
  if (n_samples < m_params.min_sample ||
      n_samples % m_params.check_frequency != 0) {
    return false;
  }
  _check(sampled_data);
  return m_results.is_complete;
}

}  // namespace Monte2
}  // namespace CASM

#endif
