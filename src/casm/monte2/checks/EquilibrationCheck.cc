#include "casm/monte2/checks/EquilibrationCheck.hh"

#include <boost/math/special_functions/erf.hpp>

namespace CASM {
namespace Monte2 {

/// \brief Check if a range of observations have equilibrated
///
/// This uses the following algorithm, based on:
///  Van de Walle and Asta, Modelling Simul. Mater. Sci. Eng. 10 (2002) 521–538.
///
/// Partition observations into three ranges:
///
///   - equilibriation stage:  [0, start1)
///   - first partition:  [start1, start2)
///   - second partition: [start2, N)
///
/// where N is observations.size(), start1 and start 2 are indices into
/// observations such 0 <= start1 < start2 <= N, the number of elements in
/// the first and second partition are the same (within 1).
///
/// The calculation is considered equilibrated at start1 if the mean of the
/// elements in the first and second partition are approximately equal to the
/// desired precsion: (std::abs(mean1 - mean2) < prec).
///
/// Additionally, the value start1 is incremented as much as needed to ensure
/// that the equilibriation stage has observations on either side of the total
/// mean.
///
/// If all observations are approximately equal, then:
/// - is_equilibrated = true
/// - N_samples_for_equilibration = 0
///
/// If the equilibration conditions are met, the result contains:
/// - is_equilibrated = true
/// - N_samples_for_equilibration = start1
///
/// If the equilibration conditions are not met, the result contains:
/// - is_equilibrated = false
/// - N_samples_for_equilibration = <undefined>
///
/// \param observations An Eigen::VectorXd of observations
/// \param prec Desired absolute precision (<X> +/- prec)
/// \param check_all If true, return results for all requested sampler
/// components. If false, return when the first component that is not
/// equilibrated is encountered.
/// \returns An IndividualEquilibrationCheckResult instance
///
IndividualEquilibrationCheckResult equilibration_check(
    Eigen::VectorXd const &observations, double prec) {
  CountType start1, start2, N;
  IndividualEquilibrationCheckResult result;
  double sum1, sum2;
  double eps =
      (observations(0) == 0.0) ? 1e-8 : std::abs(observations(0)) * 1e-8;

  N = observations.size();
  bool is_even = ((N % 2) == 0);

  // -------------------------------------------------
  // if the values are all the same to observations(0)*1e-8, set
  // m_is_equilibrated = true; m_equil_samples = 0;
  bool all_same = true;
  for (CountType i = 0; i < N; i++)
    if (std::abs(observations(i) - observations(0)) > eps) {
      all_same = false;
      break;
    }
  if (all_same) {
    result.is_equilibrated = true;
    result.N_samples_for_equilibration = 0;
    return result;
  }

  // find partitions
  start1 = 0;
  start2 = (is_even) ? N / 2 : (N / 2) + 1;

  // find sums for each partition
  sum1 = observations.head(start2).sum();
  sum2 = observations.segment(start2, N - start2).sum();

  // increment start1 (and update start2, sum1, and sum2)
  // until abs(mean1 - mean2) < prec
  while (std::abs((sum1 / (start2 - start1)) - (sum2 / (N - start2))) > prec &&
         start1 < N - 2) {
    if (is_even) {
      sum1 -= observations(start1);
      sum1 += observations(start2);
      sum2 -= observations(start2);
      start2++;
    } else {
      sum1 -= observations(start1);
    }

    start1++;
    is_even = !is_even;
  }

  // ensure that the equilibration stage
  // has values on either side of the total mean
  double mean_tot = (sum1 + sum2) / (N - start1);
  if (observations(start1) < mean_tot) {
    while (observations(start1) < mean_tot && start1 < N - 1) start1++;
  } else {
    while (observations(start1) > mean_tot && start1 < N - 1) start1++;
  }

  result.is_equilibrated = (start1 < N - 1);
  result.N_samples_for_equilibration = start1;
  return result;
}

/// \brief Check convergence of all requested properties
///
/// \param convergence_check_params Sampler components to check, with requested
///     precision and confidence
/// \param N_samples_for_equilibration Number of samples to exclude from
///     statistics because the system is out of equilibrium
/// \param samplers All samplers
///
/// \returns A ConvergenceCheckResults instance. Note that
///     N_samples_for_statistics is set to the total number of samples
///     if no convergence checks are requested (when
///     `convergence_check_params.size() == 0`), otherwise it will be equal to
///     `get_n_samples(samplers) - N_samples_for_equilibration`.
EquilibrationCheckResults equilibration_check(
    std::map<SamplerComponent, ConvergenceCheckParams> const
        &convergence_check_params,
    std::map<std::string, std::shared_ptr<Sampler>> const &samplers,
    bool check_all) {
  EquilibrationCheckResults results;

  if (!convergence_check_params.size()) {
    return results;
  }

  // will set to false if any requested sampler components are not equilibrated
  results.all_equilibrated = true;

  // check requested sampler components for equilibration
  for (auto const &p : convergence_check_params) {
    SamplerComponent const &key = p.first;
    ConvergenceCheckParams const &value = p.second;

    // find and validate sampler name && component index
    auto sampler_it = find_and_validate(key, samplers);

    // do equilibration check
    IndividualEquilibrationCheckResult current = equilibration_check(
        sampler_it->second->component(key.component_index),  // observations
        value.precision);

    // combine results
    results.N_samples_for_all_to_equilibrate =
        std::max(results.N_samples_for_all_to_equilibrate,
                 current.N_samples_for_equilibration);
    results.all_equilibrated &= current.is_equilibrated;
    results.individual_results.emplace(key, current);

    // break if possible
    if (!check_all && !results.all_equilibrated) {
      break;
    }
  }
  return results;
}

}  // namespace Monte2
}  // namespace CASM
