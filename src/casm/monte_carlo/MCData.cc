#include "casm/monte_carlo/MCData.hh"

#include "casm/external/boost.hh"

namespace CASM {

  /// \brief Check if a range of observations have equilibrated
  ///
  /// \param observations An Eigen::VectorXd of observations
  /// \param prec Desired absolute precision (<X> +/- prec)
  ///
  /// This uses the following algorithm, based on:
  ///  Van de Walle and Asta, Modelling Simul. Mater. Sci. Eng. 10 (2002) 521–538.
  ///
  /// Partition observations into three ranges:
  ///
  ///  - equilibriation stage:  [0, start1)
  ///  - first partition:  [start1, start2)
  ///  - second partition: [start2, N)
  ///
  ///    where N is observations.size(), start1 and start 2 are indices into observations
  ///    such 0 <= start1 < start2 <= N, the number of elements in the first and second
  ///    partition are the same (within 1).
  ///
  ///  We consider the calculation equilibrated at start1 if mean of elements in first and second partition
  ///  are approximately equal to the desired precsion: (std::abs(mean1 - mean2) < prec).
  ///
  ///  Additionally, we increment start1 as much as needed to ensure that the equilibriation stage has
  ///  observations on either side of the total mean.
  ///
  ///  If these conditions are met,
  ///    set: m_is_equilibrated = true; m_equil_samples = start1;
  ///
  ///  If these conditions can not be met,
  ///    set: m_is_equilibrated = false; m_equil_samples = N-1;
  ///
  MCDataEquilibration::MCDataEquilibration(const Eigen::VectorXd &observations, double prec) {

    size_type start1, start2, N;
    double sum1, sum2;
    double eps = (observations(0) == 0.0) ? 1e-8 : std::abs(observations(0)) * 1e-8;

    N = observations.size();
    bool is_even = ((N % 2) == 0);

    // -------------------------------------------------
    // if the values are all the same to observations(0)*1e-8, set m_is_equilibrated = true; m_equil_samples = 0;
    bool all_same = true;
    for(size_type i = 0; i < N; i++)
      if(std::abs(observations(i) - observations(0)) > eps) {
        all_same = false;
        break;
      }
    if(all_same) {
      m_is_equilibrated = true;
      m_equil_samples = 0;
      return;
    }

    // find partitions
    start1 = 0;
    start2 = (is_even) ? N / 2 : (N / 2) + 1;

    // find sums for each partition
    sum1 = observations.head(start2).sum();
    // mean1 = sum1 / (start2 - start1)
    sum2 = observations.segment(start2, N - start2).sum();
    // mean2 = sum2 / (N - start2)

    //std::cout << observations << std::endl;

    // increment start1 (and update start2, sum1, and sum2) until abs(mean1 - mean2) < prec
    while(std::abs((sum1 / (start2 - start1)) - (sum2 / (N - start2))) > prec && start1 < N - 2) {
      /*std::cout << "  start1: " << start1 << "  start2: " << start2 << "  N: " << N
                << "  mean1: " << (sum1 / (start2 - start1)) << "  mean2: " << (sum2 / (N - start2))
                << "  diff: " << std::abs((sum1 / (start2 - start1)) - (sum2 / (N - start2)))
                << "  prec: " << prec << std::endl;*/
      if(is_even) {
        sum1 -= observations(start1);
        sum1 += observations(start2);
        sum2 -= observations(start2);
        start2++;
      }
      else {
        sum1 -= observations(start1);
      }

      start1++;
      is_even = !is_even;
    }

    // ensure that there are values on either side of the total mean
    double mean_tot = (sum1 + sum2) / (N - start1);

    if(observations(start1) < mean_tot) {
      while(observations(start1) < mean_tot && start1 < N - 1)
        start1++;
    }
    else {
      while(observations(start1) > mean_tot && start1 < N - 1)
        start1++;
    }

    m_is_equilibrated = (start1 < N - 1);
    m_equil_samples = start1;
  }


  /// \brief Check convergence of a range of observations
  ///
  /// \param observations An Eigen::VectorXd of observations
  /// \param conf Desired confidence level
  ///
  /// We check for convergence using the algorithm of:
  ///  Van de Walle and Asta, Modelling Simul. Mater. Sci. Eng. 10 (2002) 521–538.
  ///
  /// The observations are considered converged to the desired prec and conf if:
  /// - calculated_prec <= prec,
  ///
  /// where:
  /// - calculated_prec = sqrt(Var)*z_alpha,
  /// - z_alpha = sqrt(2.0)*inv_erf(1.0-conf)
  /// - var_of_mean = (CoVar[0]/observations.size())*((1.0+rho)/(1.0-rho))
  /// - CoVar[i] = ( (1.0/(observations.size()-i))*sum_j(j=0:L-i-1, observations(j)*observations(j+1)) ) - sqr(observations.mean());
  /// - rho = pow(2.0, -1.0/i), using min i such that CoVar[i]/CoVar[0] < 0.5
  ///
  MCDataConvergence::MCDataConvergence(const Eigen::VectorXd &observations, double conf) :
    m_is_converged(false),
    m_mean(observations.mean()),
    m_squared_norm(observations.squaredNorm()) {

    // will check if Var <= criteria
    double z_alpha = sqrt(2.0) * boost::math::erf_inv(conf);

    bool found_rho;
    double rho;
    double CoVar0;
    size_type N = observations.size();

    // try to calculate variance taking into account correlations
    std::tie(found_rho, rho, CoVar0) = _calc_rho(observations);

    if(!found_rho) {
      m_is_converged = false;
      m_calculated_prec = 1.0 / 0.0;
      return;
    }


    double var_of_mean = (CoVar0 / N) * (1.0 + rho) / (1.0 - rho);

    m_calculated_prec = z_alpha * sqrt(var_of_mean);
  }

  double covariance(const Eigen::VectorXd &X, const Eigen::VectorXd &Y, double mean) {
    Eigen::VectorXd vmean = Eigen::VectorXd::Ones(X.size()) * mean;
    return (X - vmean).dot(Y - vmean) / X.size();
  }


  /// \brief Try to find rho = pow(2.0, -1.0/i), using min i such that CoVar[i]/CoVar[0] <= 0.5
  ///
  /// \returns std::tuple<bool, double> : (found_rho?, rho)
  ///
  std::tuple<bool, double, double> MCDataConvergence::_calc_rho(const Eigen::VectorXd &observations) {

    size_type N = observations.size();
    double CoVar0 = m_squared_norm / N - m_mean * m_mean;

    // if there is essentially no variation, return rho(l==1)
    if(std::abs(CoVar0 / m_mean) < 1e-8 || CoVar0 == 0.0) {
      CoVar0 = 0.0;
      return std::make_tuple(true, pow(2.0, -1.0 / 1), CoVar0);
    }

    // simple incremental search for now, could try bracketing / bisection
    for(size_type i = 1; i < N; ++i) {

      size_type range_size = N - i;

      double cov = covariance(observations.segment(0, range_size), observations.segment(i, range_size), m_mean);

      if(std::abs(cov / CoVar0) <= 0.5) {
        return std::make_tuple(true, pow(2.0, (-1.0 / i)), CoVar0);
      }
    }

    // if could not find:
    return std::make_tuple(false, 0.0, CoVar0);
  }

}
