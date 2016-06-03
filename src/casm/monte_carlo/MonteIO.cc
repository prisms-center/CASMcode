#include "casm/monte_carlo/MonteIO.hh"

namespace CASM {

  /// \brief Print mean property values: <prop_name>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloMeanFormatter(std::string prop_name) {

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      return mc->samplers().find(prop_name)->second->mean(mc->is_equilibrated().second);
    };

    auto validator = [ = ](const ConstMonteCarloPtr & mc) {
      return mc->is_equilibrated().first;
    };

    std::string header = std::string("<") + prop_name + ">";

    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator, validator);
  }

  /// \brief Print calculated precision of property values: prec(<prop_name>)
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloPrecFormatter(std::string prop_name) {

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      return mc->samplers().find(prop_name)->second->calculated_precision(mc->is_equilibrated().second);
    };

    auto validator = [ = ](const ConstMonteCarloPtr & mc) {
      return mc->is_equilibrated().first;
    };

    std::string header = std::string("prec(<") + prop_name + ">)";

    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator, validator);
  }

  double CovEvaluator::operator()(const ConstMonteCarloPtr &mc) {
    auto equil = mc->is_equilibrated();

    // cov = <X*Y> - <X>*<Y>
    const MonteSampler &sampler1 = *(mc->samplers().find(prop_name1)->second);
    const Eigen::VectorXd &obs1 = sampler1.data().observations();
    const Eigen::VectorXd &X = obs1.segment(equil.second, obs1.size() - equil.second);
    double k1 = X(0);

    const MonteSampler &sampler2 = *(mc->samplers().find(prop_name2)->second);
    const Eigen::VectorXd &obs2 = sampler2.data().observations();
    const Eigen::VectorXd &Y = obs2.segment(equil.second, obs2.size() - equil.second);
    double k2 = Y(0);

    double Xsum = 0.0;
    double Ysum = 0.0;
    double XYsum = 0.0;
    Index N = X.size();

    for(Index i = 0; i < N; ++i) {
      Xsum += X(i);
      Ysum += Y(i);
      XYsum += X(i) * Y(i);
    }

    //return X.cwiseProduct(Y).mean() - X.mean()*Y.mean();
    return (XYsum - Xsum * Ysum / N) / N;
  }

  /// \brief Print covariance: cov(prop_name1, prop_name2)
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloCovFormatter(std::string prop_name1, std::string prop_name2) {

    auto evaluator = CovEvaluator(prop_name1, prop_name2);

    auto validator = [ = ](const ConstMonteCarloPtr & mc) {
      return mc->is_equilibrated().first;
    };

    std::string header = std::string("cov(") + prop_name1 + "," + prop_name2 + ")";

    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator, validator);
  }

  /// \brief Print if equilibrated (not counting explicitly requested equilibration)
  GenericDatumFormatter<bool, ConstMonteCarloPtr> MonteCarloIsEquilibratedFormatter() {

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc)->bool {
      return mc->is_equilibrated().first;
    };

    return GenericDatumFormatter<bool, ConstMonteCarloPtr>("is_equilibrated", "is_equilibrated", evaluator);
  }

  /// \brief Print if converged
  GenericDatumFormatter<bool, ConstMonteCarloPtr> MonteCarloIsConvergedFormatter() {

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc)->bool {
      return mc->is_converged();
    };

    return GenericDatumFormatter<bool, ConstMonteCarloPtr>("is_converged", "is_converged", evaluator);
  }

  /// \brief Print number of samples used for equilibration (not counting explicitly requested equilibration)
  GenericDatumFormatter<MonteSampler::size_type, ConstMonteCarloPtr> MonteCarloNEquilSamplesFormatter() {

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc)->unsigned long {
      return mc->is_equilibrated().second;
    };

    return GenericDatumFormatter<MonteSampler::size_type, ConstMonteCarloPtr>("N_equil_samples", "N_equil_samples", evaluator);
  }

  /// \brief Print number of samples used in calculating means
  GenericDatumFormatter<MonteSampler::size_type, ConstMonteCarloPtr> MonteCarloNAvgSamplesFormatter() {

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc)->unsigned long {
      return mc->sample_times().size() - mc->is_equilibrated().second;
    };

    return GenericDatumFormatter<MonteSampler::size_type, ConstMonteCarloPtr>("N_avg_samples", "N_avg_samples", evaluator);
  }

  /// \brief Print Pass number of observation
  GenericDatumFormatter<MonteCounter::size_type, std::pair<ConstMonteCarloPtr, Index> > MonteCarloPassFormatter() {

    auto evaluator = [ = ](const std::pair<ConstMonteCarloPtr, Index> &obs) {
      return obs.first->sample_times()[obs.second].first;
    };

    return GenericDatumFormatter<MonteCounter::size_type, std::pair<ConstMonteCarloPtr, Index> >("Pass", "Pass", evaluator);
  }

  /// \brief Print Step number of observation
  GenericDatumFormatter<MonteCounter::size_type, std::pair<ConstMonteCarloPtr, Index> > MonteCarloStepFormatter() {

    auto evaluator = [ = ](const std::pair<ConstMonteCarloPtr, Index> &obs) {
      return obs.first->sample_times()[obs.second].second;
    };

    return GenericDatumFormatter<MonteCounter::size_type, std::pair<ConstMonteCarloPtr, Index> >("Step", "Step", evaluator);
  }

  /// \brief Print value of observation
  GenericDatumFormatter<double, std::pair<ConstMonteCarloPtr, Index> > MonteCarloObservationFormatter(std::string prop_name) {

    auto evaluator = [ = ](const std::pair<ConstMonteCarloPtr, Index> &obs) {
      const MonteSampler &sampler = *(obs.first->samplers().find(prop_name)->second);
      return sampler.data().observations()(obs.second);
    };

    return GenericDatumFormatter<double, std::pair<ConstMonteCarloPtr, Index> >(prop_name, prop_name, evaluator);
  }

  /// \brief Print value of a particular occupation variable
  GenericDatumFormatter<int, std::pair<ConstMonteCarloPtr, Index> > MonteCarloOccFormatter(Index occ_index) {

    auto evaluator = [ = ](const std::pair<ConstMonteCarloPtr, Index> &site)->int {
      return site.first->trajectory()[site.second].occ(occ_index);
    };

    std::string header = std::string("occ(") + std::to_string(occ_index) + ")";

    return GenericDatumFormatter<int, std::pair<ConstMonteCarloPtr, Index> >(header, header, evaluator);
  }


}

