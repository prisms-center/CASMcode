#include "casm/monte2/checks/CompletionCheck.hh"

#include <boost/math/special_functions/erf.hpp>

namespace CASM {
namespace Monte2 {

CompletionCheck::CompletionCheck(CompletionCheckParams params)
    : m_params(params) {}

/// \brief Check for equilibration and convergence, then set m_results
void CompletionCheck::_check(
    std::map<std::string, std::shared_ptr<Sampler>> const &samplers,
    std::optional<CountType> count, std::optional<TimeType> time,
    CountType n_samples) {
  m_results = CompletionCheckResults();

  // if minimums not met -> continue
  if (!all_minimums_met(m_params.cutoff_params, count, time, n_samples)) {
    m_results.is_complete = false;
    return;
  }

  // if auto convergence mode:
  if (m_params.convergence_check_params.size()) {
    // check for equilibration
    bool check_all = false;
    m_results.equilibration_check_results = equilibration_check(
        m_params.convergence_check_params, samplers, check_all);

    // if all requested to converge are equilibrated, then check convergence
    if (m_results.equilibration_check_results.all_equilibrated) {
      m_results.convergence_check_results =
          convergence_check(m_params.convergence_check_params,
                            m_results.equilibration_check_results
                                .N_samples_for_all_to_equilibrate,
                            samplers);
    }

    // if all requested to converge are converged, then complete
    if (m_results.convergence_check_results.all_converged) {
      m_results.is_complete = true;
      return;
    }
  }

  // if any maximum met, stop
  if (any_maximum_met(m_params.cutoff_params, count, time, n_samples)) {
    m_results.is_complete = true;
  }
}

}  // namespace Monte2
}  // namespace CASM
