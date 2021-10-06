#ifndef CASM_monte2_CompletionCheck
#define CASM_monte2_CompletionCheck

#include <optional>

#include "casm/monte2/checks/ConvergenceCheck.hh"
#include "casm/monte2/checks/CutoffCheck.hh"
#include "casm/monte2/checks/EquilibrationCheck.hh"
#include "casm/monte2/definitions.hh"
#include "casm/monte2/sampling/SampledData.hh"
#include "casm/monte2/sampling/Sampler.hh"

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
  CountType check_begin = 0;

  /// \brief How often to check for completion
  ///
  /// Check for completion performed if:
  /// - n_samples % check_frequency == 0 && n_samples >= check_begin
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

  bool is_complete(
      std::map<std::string, std::shared_ptr<Sampler>> const &samplers);

  bool is_complete(
      std::map<std::string, std::shared_ptr<Sampler>> const &samplers,
      CountType count);

  bool is_complete(
      std::map<std::string, std::shared_ptr<Sampler>> const &samplers,
      TimeType time);

  bool is_complete(
      std::map<std::string, std::shared_ptr<Sampler>> const &samplers,
      CountType count, TimeType time);

  bool is_complete(SampledData const &sampled_data);

  CompletionCheckResults const &results() const { return m_results; }

 private:
  bool _is_complete(
      std::map<std::string, std::shared_ptr<Sampler>> const &samplers,
      std::optional<CountType> count, std::optional<TimeType> time);

  void _check(std::map<std::string, std::shared_ptr<Sampler>> const &samplers,
              std::optional<CountType> count, std::optional<TimeType> time,
              CountType n_samples);

  CompletionCheckParams m_params;

  CompletionCheckResults m_results;
};

// --- Inline definitions ---

inline bool CompletionCheck::is_complete(
    std::map<std::string, std::shared_ptr<Sampler>> const &samplers) {
  return _is_complete(samplers, std::nullopt, std::nullopt);
}

inline bool CompletionCheck::is_complete(
    std::map<std::string, std::shared_ptr<Sampler>> const &samplers,
    CountType count) {
  return _is_complete(samplers, count, std::nullopt);
}

inline bool CompletionCheck::is_complete(
    std::map<std::string, std::shared_ptr<Sampler>> const &samplers,
    TimeType time) {
  return _is_complete(samplers, std::nullopt, time);
}

inline bool CompletionCheck::is_complete(
    std::map<std::string, std::shared_ptr<Sampler>> const &samplers,
    CountType count, TimeType time) {
  return _is_complete(samplers, count, time);
}

inline bool CompletionCheck::is_complete(SampledData const &sampled_data) {
  return _is_complete(sampled_data.samplers, get_count(sampled_data),
                      get_time(sampled_data));
}

inline bool CompletionCheck::_is_complete(
    std::map<std::string, std::shared_ptr<Sampler>> const &samplers,
    std::optional<CountType> count, std::optional<TimeType> time) {
  CountType n_samples = get_n_samples(samplers);
  if (n_samples < m_params.check_begin ||
      n_samples % m_params.check_frequency != 0) {
    return false;
  }
  _check(samplers, count, time, n_samples);
  return m_results.is_complete;
}

}  // namespace Monte2
}  // namespace CASM

#endif
