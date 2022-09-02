#include "casm/monte_carlo/CorrMatchingPotential.hh"

#include <iostream>

#include "casm/clex/random_alloy_correlations_impl.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"

namespace CASM {
namespace Monte {

double corr_matching_potential(Eigen::VectorXd const &corr,
                               CorrMatchingParams const &params) {
  double Epot = 0;
  bool counting_n_exact = true;
  Index n_exact = 0;

  for (auto const &target : params.targets) {
    if (target.index < 0 || target.index >= corr.size()) {
      throw std::runtime_error(
          "Error calculating correlations matching potential: target index out "
          "of range");
    }
    double value = corr(target.index);
    if (counting_n_exact) {
      if (CASM::almost_equal(value, target.value, params.tol)) {
        ++n_exact;
      } else {
        counting_n_exact = false;
      }
    }
    Epot += target.weight * std::abs(value - target.value);
  }
  Epot -= params.exact_matching_weight * n_exact;
  return Epot;
}

Eigen::VectorXd make_corr_matching_error(Eigen::VectorXd const &corr,
                                         CorrMatchingParams const &params) {
  Eigen::VectorXd corr_matching_error = Eigen::VectorXd::Zero(corr.size());

  for (auto const &target : params.targets) {
    if (target.index < 0 || target.index >= corr.size()) {
      throw std::runtime_error(
          "Error calculating correlations matching potential: target index out "
          "of range");
    }
    corr_matching_error(target.index) = corr(target.index) - target.value;
  }
  return corr_matching_error;
}

double delta_corr_matching_potential(Eigen::VectorXd const &corr,
                                     Eigen::VectorXd const &delta_corr,
                                     CorrMatchingParams const &params) {
  if (corr.size() != delta_corr.size()) {
    throw std::runtime_error(
        "Error calculating correlations matching potential delta: corr and "
        "delta_corr size mismatch");
  }
  double dEpot = 0;
  bool counting_n_exact_1 = true;
  Index n_exact_1 = 0;
  bool counting_n_exact_2 = true;
  Index n_exact_2 = 0;
  for (auto const &target : params.targets) {
    if (target.index < 0 || target.index >= corr.size()) {
      throw std::runtime_error(
          "Error calculating correlations matching potential delta: target "
          "index out of range");
    }
    double value = corr(target.index);
    double dvalue = delta_corr(target.index);
    if (counting_n_exact_1) {
      if (CASM::almost_equal(value, target.value, params.tol)) {
        ++n_exact_1;
      } else {
        counting_n_exact_1 = false;
      }
    }
    if (counting_n_exact_2) {
      if (CASM::almost_equal(value + dvalue, target.value, params.tol)) {
        ++n_exact_2;
      } else {
        counting_n_exact_2 = false;
      }
    }
    dEpot += target.weight * (std::abs(value + dvalue - target.value) -
                              std::abs(value - target.value));
  }
  dEpot -= params.exact_matching_weight * (n_exact_2 - n_exact_1);
  return dEpot;
}

void incr(std::optional<CorrMatchingParams> &lhs,
          std::optional<CorrMatchingParams> const &rhs) {
  if (lhs.has_value() != rhs.has_value()) {
    throw std::runtime_error(
        "Error incrementing correlation-matching params: existence mismatch");
  }
  if (!lhs.has_value()) {
    return;
  }
  if (lhs->targets.size() != rhs->targets.size()) {
    throw std::runtime_error(
        "Error incrementing correlation-matching params: size mismatch");
  }
  for (Index i = 0; i < lhs->targets.size(); ++i) {
    if (lhs->targets[i].index != rhs->targets[i].index) {
      throw std::runtime_error(
          "Error incrementing correlation-matching params: index mismatch");
    }
    lhs->targets[i].value += rhs->targets[i].value;
    lhs->targets[i].weight += rhs->targets[i].weight;
  }
  lhs->exact_matching_weight += rhs->exact_matching_weight;
}

void decr(std::optional<CorrMatchingParams> &lhs,
          std::optional<CorrMatchingParams> const &rhs) {
  if (lhs.has_value() != rhs.has_value()) {
    throw std::runtime_error(
        "Error decrementing correlation-matching params: existence mismatch");
  }
  if (!lhs.has_value()) {
    return;
  }
  if (lhs->targets.size() != rhs->targets.size()) {
    throw std::runtime_error(
        "Error decrementing correlation-matching params: size mismatch");
  }
  for (Index i = 0; i < lhs->targets.size(); ++i) {
    if (lhs->targets[i].index != rhs->targets[i].index) {
      throw std::runtime_error(
          "Error decrementing correlation-matching params: index mismatch");
    }
    lhs->targets[i].value -= rhs->targets[i].value;
    lhs->targets[i].weight -= rhs->targets[i].weight;
  }
  lhs->exact_matching_weight -= rhs->exact_matching_weight;
}

bool almost_equal(std::optional<CorrMatchingParams> const &lhs,
                  std::optional<CorrMatchingParams> const &rhs, double tol) {
  if (lhs.has_value() != rhs.has_value()) {
    throw std::runtime_error(
        "Error checking correlation-matching params equivalence: existence "
        "mismatch");
  }
  if (!lhs.has_value()) {
    return true;
  }
  if (lhs->targets.size() != rhs->targets.size()) {
    throw std::runtime_error(
        "Error checking correlation-matching params equivalence: size "
        "mismatch");
  }
  for (Index i = 0; i < lhs->targets.size(); ++i) {
    if (lhs->targets[i].index != rhs->targets[i].index) {
      throw std::runtime_error(
          "Error checking correlation-matching params equivalence: index "
          "mismatch");
    }
    if (!CASM::almost_equal(lhs->targets[i].value, rhs->targets[i].value,
                            tol)) {
      return false;
    }
    if (!CASM::almost_equal(lhs->targets[i].weight, rhs->targets[i].weight,
                            tol)) {
      return false;
    }
  }
  return CASM::almost_equal(lhs->exact_matching_weight,
                            rhs->exact_matching_weight, tol);
}

/// \brief If lhs and rhs have values, update max_division if round((*lhs)(i) /
/// (*rhs)(i)) > max_division for any i, and ignoring divide by 0
void find_max_division(int &max_division,
                       std::optional<CorrMatchingParams> const &lhs,
                       std::optional<CorrMatchingParams> const &rhs) {
  if (lhs.has_value() != rhs.has_value()) {
    throw std::runtime_error(
        "Error finding correlation-matching params division: existence "
        "mismatch");
  }
  if (!lhs.has_value()) {
    return;
  }
  if (lhs->targets.size() != rhs->targets.size()) {
    throw std::runtime_error(
        "Error finding correlation-matching params division: size "
        "mismatch");
  }

  auto _check = [&](double lhs_value, double rhs_value) {
    if (CASM::almost_zero(rhs_value)) {
      return;
    }
    int temp_division = round(lhs_value / rhs_value);
    if (temp_division > max_division) {
      max_division = temp_division;
    }
  };

  for (Index i = 0; i < lhs->targets.size(); ++i) {
    if (lhs->targets[i].index != rhs->targets[i].index) {
      throw std::runtime_error(
          "Error checking correlation-matching params equivalence: index "
          "mismatch");
    }
    _check(lhs->targets[i].value, rhs->targets[i].value);
    _check(lhs->targets[i].weight, rhs->targets[i].weight);
  }
  _check(lhs->exact_matching_weight, rhs->exact_matching_weight);
}

void print_param(std::ostream &sout, std::string name,
                 std::optional<CorrMatchingParams> const &params) {
  if (!params.has_value()) {
    return;
  }
  sout << name << ".exact_matching_weight: " << params->exact_matching_weight
       << "\n";
  for (auto const &target : params->targets) {
    sout << name << ".target: {index=" << target.index
         << ", value=" << target.value << ", weight=" << target.weight << "}\n";
  }
}

RandomAlloyCorrMatchingParams::RandomAlloyCorrMatchingParams() {}

RandomAlloyCorrMatchingParams::RandomAlloyCorrMatchingParams(
    std::vector<Eigen::VectorXd> const &_sublattice_prob,
    std::shared_ptr<RandomAlloyCorrCalculator> const &_random_alloy_corr_f,
    double _exact_matching_weight, double _tol)
    : CorrMatchingParams(),
      sublattice_prob(_sublattice_prob),
      random_alloy_corr_f(_random_alloy_corr_f) {
  this->exact_matching_weight = _exact_matching_weight;
  this->tol = _tol;
  update_targets();
}

void RandomAlloyCorrMatchingParams::update_targets() {
  Eigen::VectorXd random_alloy_corr = (*random_alloy_corr_f)(sublattice_prob);
  this->targets.clear();
  for (int i = 0; i < random_alloy_corr.size(); ++i) {
    this->targets.push_back(CorrMatchingTarget(i, random_alloy_corr(i), 1.0));
  }
}

void incr(std::optional<RandomAlloyCorrMatchingParams> &lhs,
          std::optional<RandomAlloyCorrMatchingParams> const &rhs) {
  if (lhs.has_value() != rhs.has_value()) {
    throw std::runtime_error(
        "Error incrementing random-alloy-correlation-matching params: "
        "existence mismatch");
  }
  if (!lhs.has_value()) {
    return;
  }
  std::runtime_error e(
      "Error incrementing random-alloy-correlation-matching params: size "
      "mismatch");

  if (lhs->sublattice_prob.size() != rhs->sublattice_prob.size()) {
    throw e;
  }
  for (Index i = 0; i < lhs->sublattice_prob.size(); ++i) {
    if (lhs->sublattice_prob[i].size() != rhs->sublattice_prob[i].size()) {
      throw e;
    }
    for (Index j = 0; j < lhs->sublattice_prob[i].size(); ++j) {
      lhs->sublattice_prob[i](j) += rhs->sublattice_prob[i](j);
    }
  }
  lhs->exact_matching_weight += rhs->exact_matching_weight;
  lhs->update_targets();
}

void decr(std::optional<RandomAlloyCorrMatchingParams> &lhs,
          std::optional<RandomAlloyCorrMatchingParams> const &rhs) {
  if (lhs.has_value() != rhs.has_value()) {
    throw std::runtime_error(
        "Error decrementing random-alloy-correlation-matching params: "
        "existence mismatch");
  }
  if (!lhs.has_value()) {
    return;
  }
  std::runtime_error e(
      "Error decrementing random-alloy-correlation-matching params: size "
      "mismatch");

  if (lhs->sublattice_prob.size() != rhs->sublattice_prob.size()) {
    throw e;
  }
  for (Index i = 0; i < lhs->sublattice_prob.size(); ++i) {
    if (lhs->sublattice_prob[i].size() != rhs->sublattice_prob[i].size()) {
      throw e;
    }
    for (Index j = 0; j < lhs->sublattice_prob[i].size(); ++j) {
      lhs->sublattice_prob[i](j) -= rhs->sublattice_prob[i](j);
    }
  }
  lhs->exact_matching_weight -= rhs->exact_matching_weight;
  lhs->update_targets();
}

bool almost_equal(std::optional<RandomAlloyCorrMatchingParams> const &lhs,
                  std::optional<RandomAlloyCorrMatchingParams> const &rhs,
                  double tol) {
  if (lhs.has_value() != rhs.has_value()) {
    throw std::runtime_error(
        "Error checking random-alloy-correlation-matching params equivalence: "
        "existence "
        "mismatch");
  }
  if (!lhs.has_value()) {
    return true;
  }
  std::runtime_error e(
      "Error checking random-alloy-correlation-matching params equivalence: "
      "size "
      "mismatch");

  if (lhs->sublattice_prob.size() != rhs->sublattice_prob.size()) {
    throw e;
  }

  for (Index i = 0; i < lhs->sublattice_prob.size(); ++i) {
    if (lhs->sublattice_prob[i].size() != rhs->sublattice_prob[i].size()) {
      throw e;
    }
    for (Index j = 0; j < lhs->sublattice_prob[i].size(); ++j) {
      if (!CASM::almost_equal(lhs->sublattice_prob[i](j),
                              rhs->sublattice_prob[i](j), tol)) {
        return false;
      }
    }
  }
  return CASM::almost_equal(lhs->exact_matching_weight,
                            rhs->exact_matching_weight, tol);
}

/// \brief If lhs and rhs have values, update max_division if round((*lhs)(i) /
/// (*rhs)(i)) > max_division for any i, and ignoring divide by 0
void find_max_division(
    int &max_division, std::optional<RandomAlloyCorrMatchingParams> const &lhs,
    std::optional<RandomAlloyCorrMatchingParams> const &rhs) {
  if (lhs.has_value() != rhs.has_value()) {
    throw std::runtime_error(
        "Error finding random-alloy-correlation-matching params divisions: "
        "existence "
        "mismatch");
  }
  if (!lhs.has_value()) {
    return;
  }
  std::runtime_error e(
      "Error finding random-alloy-correlation-matching params divisions: size "
      "mismatch");

  if (lhs->sublattice_prob.size() != rhs->sublattice_prob.size()) {
    throw e;
  }

  auto _check = [&](double lhs_value, double rhs_value) {
    if (CASM::almost_zero(rhs_value)) {
      return;
    }
    int temp_division = round(lhs_value / rhs_value);
    if (temp_division > max_division) {
      max_division = temp_division;
    }
  };

  for (Index i = 0; i < lhs->sublattice_prob.size(); ++i) {
    if (lhs->sublattice_prob[i].size() != rhs->sublattice_prob[i].size()) {
      throw e;
    }
    for (Index j = 0; j < lhs->sublattice_prob[i].size(); ++j) {
      _check(lhs->sublattice_prob[i](j), rhs->sublattice_prob[i](j));
    }
  }
  _check(lhs->exact_matching_weight, rhs->exact_matching_weight);
}

void print_param(std::ostream &sout, std::string name,
                 std::optional<RandomAlloyCorrMatchingParams> const &params) {
  if (!params.has_value()) {
    return;
  }
  sout << name << ".sublattice_prob: " << std::endl;
  Index b = 0;
  for (auto const &prob : params->sublattice_prob) {
    sout << "  " << b << ": " << prob.transpose() << std::endl;
    ++b;
  }
  sout << name << ".exact_matching_weight: " << params->exact_matching_weight
       << "\n";
}

}  // namespace Monte
}  // namespace CASM
