#ifndef CASM_monte2_ConvergenceCheck
#define CASM_monte2_ConvergenceCheck

#include "casm/monte2/Definitions.hh"
#include "casm/monte2/Sampling.hh"

namespace CASM {
namespace Monte2 {

/// \brief Check sampler of given name for equilibration and convergence
struct ConvergenceCheckParams {
  /// \brief Requested precision
  double precision = 0.0;

  /// \brief Requested confidence
  double confidence = 0.0;
};

/// \brief Convergence check results data structure (individual component)
struct IndividualConvergenceCheckResult {
  /// \brief True if mean (<X>) is converged to requested precision
  bool is_converged;

  /// \brief Mean of property (<X>)
  double mean;

  /// \brief Squared norm of property (<X*X>)
  double squared_norm;

  /// \brief Calculated absolute precision in <X>
  ///
  /// Notes:
  /// - See `convergence_check` function for calculation details
  double calculated_precision;
};

/// \brief Check convergence of a range of observations
IndividualConvergenceCheckResult convergence_check(
    Eigen::VectorXd const &observations, double precision, double confidence);

/// \brief Convergence check results data structure (all requested components)
struct ConvergenceCheckResults {
  /// \brief True if all required properties are converged to required precision
  ///
  /// Notes:
  /// - True if completion check finds all required properties are converged to
  /// the requested precision
  /// - False otherwise, including if no convergence checks were requested
  bool all_converged = false;

  /// \brief How many samples were used to get statistics
  ///
  /// Notes:
  /// - Set to the total number of samples if no convergence checks were
  /// requested
  CountType N_samples_for_statistics = 0;

  /// \brief Results from checking convergence criteria
  std::map<SamplerComponent, IndividualConvergenceCheckResult>
      individual_results;
};

/// \brief Check convergence of all requested properties
ConvergenceCheckResults convergence_check(
    std::map<SamplerComponent, ConvergenceCheckParams> const
        &convergence_check_params,
    CountType N_samples_for_equilibration,
    std::map<std::string, Sampler> const &samplers);

}  // namespace Monte2
}  // namespace CASM

#endif
