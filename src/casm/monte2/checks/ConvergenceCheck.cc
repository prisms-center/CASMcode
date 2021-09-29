#include "casm/monte2/checks/ConvergenceCheck.hh"

#include <boost/math/special_functions/erf.hpp>

namespace CASM {
namespace Monte2 {

namespace {  // anonymous

double _covariance(const Eigen::VectorXd &X, const Eigen::VectorXd &Y,
                   double mean) {
  Eigen::VectorXd vmean = Eigen::VectorXd::Ones(X.size()) * mean;
  return (X - vmean).dot(Y - vmean) / X.size();
}

/// \brief Try to find rho = pow(2.0, -1.0/i), using min i such that
/// CoVar[i]/CoVar[0] <= 0.5
///
/// \returns std::tuple<bool, double, double> : (found_rho?, rho, CoVar[0])
///
std::tuple<bool, double, double> _calc_rho(Eigen::VectorXd const &observations,
                                           double mean, double squared_norm) {
  CountType N = observations.size();
  double CoVar0 = squared_norm / N - mean * mean;

  // if there is essentially no variation, return rho(l==1)
  if (std::abs(CoVar0 / mean) < 1e-8 || CoVar0 == 0.0) {
    CoVar0 = 0.0;
    return std::make_tuple(true, pow(2.0, -1.0 / 1), CoVar0);
  }

  // simple incremental search for now, could try bracketing / bisection
  for (CountType i = 1; i < N; ++i) {
    CountType range_size = N - i;

    double cov = _covariance(observations.segment(0, range_size),
                             observations.segment(i, range_size), mean);

    if (std::abs(cov / CoVar0) <= 0.5) {
      return std::make_tuple(true, pow(2.0, (-1.0 / i)), CoVar0);
    }
  }

  // if could not find:
  return std::make_tuple(false, 0.0, CoVar0);
}

}  // namespace

/// \brief Check convergence of a range of observations
///
/// Convergence is checked using the algorithm of:
///  Van de Walle and Asta, Modelling Simul. Mater. Sci. Eng. 10 (2002) 521–538.
///
/// The observations are considered converged to the desired precision at a
/// particular confidence level if:
///
///     calculated_precision <= precision,
///
/// where:
/// - calculated_precision = z_alpha*sqrt(var_of_mean),
/// - z_alpha = sqrt(2.0)*inv_erf(1.0-conf)
/// - var_of_mean = (CoVar[0]/observations.size())*((1.0+rho)/(1.0-rho))
/// - CoVar[i] = ( (1.0/(observations.size()-i))*sum_j(j=0:L-i-1,
/// observations(j)*observations(j+1)) ) - sqr(observations.mean());
/// - rho = pow(2.0, -1.0/i), using min i such that CoVar[i]/CoVar[0] < 0.5
///
/// \param observations An Eigen::VectorXd of observations. Should only include
///     those outside of the equilibration range.
/// \param precision Desired precision level (absolute <X> +/- precision)
/// \param confidence Desired confidence level
///
/// \returns An IndividualConvergenceCheckResult instance
///
IndividualConvergenceCheckResult convergence_check(
    Eigen::VectorXd const &observations, double precision, double confidence) {
  IndividualConvergenceCheckResult result;
  result.is_converged = false;
  result.mean = observations.mean();
  result.squared_norm = observations.squaredNorm();

  // try to calculate variance taking into account correlations
  bool found_rho;
  double rho;
  double CoVar0;
  std::tie(found_rho, rho, CoVar0) =
      _calc_rho(observations, result.mean, result.squared_norm);

  if (!found_rho) {
    result.is_converged = false;
    result.calculated_precision = 1.0 / 0.0;
    return result;
  }

  // calculated precision:
  CountType N = observations.size();
  double z_alpha = sqrt(2.0) * boost::math::erf_inv(confidence);
  double var_of_mean = (CoVar0 / N) * (1.0 + rho) / (1.0 - rho);
  result.calculated_precision = z_alpha * sqrt(var_of_mean);

  // compare to requested precision
  result.is_converged = result.calculated_precision < precision;
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
ConvergenceCheckResults convergence_check(
    std::map<SamplerComponent, ConvergenceCheckParams> const
        &convergence_check_params,
    CountType N_samples_for_equilibration,
    std::map<std::string, std::shared_ptr<Sampler>> const &samplers) {
  ConvergenceCheckResults results;
  CountType N_samples = get_n_samples(samplers);

  // if nothing to check, return
  if (!convergence_check_params.size()) {
    results.N_samples_for_statistics = N_samples;
    return results;
  }

  // if no samples after equilibration, return
  if (N_samples_for_equilibration >= N_samples) {
    return results;
  }

  // set number of samples used for statistics
  results.N_samples_for_statistics = N_samples - N_samples_for_equilibration;

  // will set to false if any requested sampler components are not converged
  results.all_converged = true;

  // check requested sampler components for equilibration
  for (auto const &p : convergence_check_params) {
    SamplerComponent const &key = p.first;
    ConvergenceCheckParams const &value = p.second;

    // find and validate sampler name && component index
    auto sampler_it = find_and_validate(key, samplers);

    // do convergence check
    IndividualConvergenceCheckResult current = convergence_check(
        sampler_it->second->component(key.component_index),  // observations
        value.precision, value.confidence);

    // combine results
    results.all_converged &= current.is_converged;
    results.individual_results.emplace(key, current);
  }
  return results;
}

}  // namespace Monte2
}  // namespace CASM
