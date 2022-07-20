#ifndef CASM_Monte_CorrMatchingPotential
#define CASM_Monte_CorrMatchingPotential

#include <optional>
#include <vector>

#include "casm/clex/random_alloy_correlations.hh"
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {
namespace Monte {

struct CorrMatchingTarget {
  CorrMatchingTarget() : index(0), value(0.0), weight(1.0) {}

  CorrMatchingTarget(Index _index, double _value, double _weight)
      : index(_index), value(_value), weight(_weight) {}

  /// \brief Correlation index
  Index index;

  /// \brief Correlation target value
  double value;

  /// \brief Weight given to difference from target
  double weight;
};

/// \brief Parameters for a correlation-matching potential
///
/// Implements:
///
///     Epot = -w_{exact}*N_{exact} +
///         \sum_i v_i * | \Gamma_{j_i} - \Gamma^{target}_{j_i}) |,
///
/// where:
/// - N_{exact} is that maximum value such that
///   \Gamma_{j_i} - \Gamma^{target}_{j_i}) ~ 0 for all i < N_{exact}
/// - exact_matching_weight = w_{exact},
/// - targets[i].index = w_i
/// - targets[i].value = \Gamma^{target}_{j_i}
/// - targets[i].weight = v_i
struct CorrMatchingParams {
  CorrMatchingParams() : exact_matching_weight(0.0), tol(CASM::TOL) {}

  CorrMatchingParams(double _exact_matching_weight, double _tol,
                     std::vector<CorrMatchingTarget> const &_targets)
      : exact_matching_weight(_exact_matching_weight),
        targets(_targets),
        tol(_tol) {}

  /// \brief Bias given for leading exactly matching correlations
  double exact_matching_weight;

  /// \brief Correlation matching targets
  std::vector<CorrMatchingTarget> targets;

  /// \brief Tolerance used to check for exactly matching correlations
  double tol;
};

double corr_matching_potential(Eigen::VectorXd const &corr,
                               CorrMatchingParams const &params);

Eigen::VectorXd make_corr_matching_error(Eigen::VectorXd const &corr,
                                         CorrMatchingParams const &params);

double delta_corr_matching_potential(Eigen::VectorXd const &corr,
                                     Eigen::VectorXd const &delta_corr,
                                     CorrMatchingParams const &params);

void incr(std::optional<CorrMatchingParams> &lhs,
          std::optional<CorrMatchingParams> const &rhs);

void decr(std::optional<CorrMatchingParams> &lhs,
          std::optional<CorrMatchingParams> const &rhs);

bool almost_equal(std::optional<CorrMatchingParams> const &lhs,
                  std::optional<CorrMatchingParams> const &rhs, double tol);

void find_max_division(int &max_division,
                       std::optional<CorrMatchingParams> const &lhs,
                       std::optional<CorrMatchingParams> const &rhs);

void print_param(std::ostream &sout, std::string name,
                 std::optional<CorrMatchingParams> const &params);

struct RandomAlloyCorrMatchingParams : public CorrMatchingParams {
  RandomAlloyCorrMatchingParams();

  RandomAlloyCorrMatchingParams(
      std::vector<Eigen::VectorXd> const &_sublattice_prob,
      std::shared_ptr<RandomAlloyCorrCalculator> const &_random_alloy_corr_f,
      double _exact_matching_weight, double _tol);

  void update_targets();

  std::vector<Eigen::VectorXd> sublattice_prob;

  std::shared_ptr<RandomAlloyCorrCalculator> random_alloy_corr_f;
};

void incr(std::optional<RandomAlloyCorrMatchingParams> &lhs,
          std::optional<RandomAlloyCorrMatchingParams> const &rhs);

void decr(std::optional<RandomAlloyCorrMatchingParams> &lhs,
          std::optional<RandomAlloyCorrMatchingParams> const &rhs);

bool almost_equal(std::optional<RandomAlloyCorrMatchingParams> const &lhs,
                  std::optional<RandomAlloyCorrMatchingParams> const &rhs,
                  double tol);

void find_max_division(int &max_division,
                       std::optional<RandomAlloyCorrMatchingParams> const &lhs,
                       std::optional<RandomAlloyCorrMatchingParams> const &rhs);

void print_param(std::ostream &sout, std::string name,
                 std::optional<RandomAlloyCorrMatchingParams> const &params);

}  // namespace Monte
}  // namespace CASM

#endif
