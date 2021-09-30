#ifndef CASM_monte2_StateSampler
#define CASM_monte2_StateSampler

#include "casm/monte2/sampling/Sampler.hh"
#include "casm/monte2/state/State.hh"

namespace CASM {
namespace Monte2 {

template <typename _ConfigType>
struct StateSamplingFunction {
  typedef _ConfigType ConfigType;

  StateSamplingFunction(
      std::string _name, std::string _description,
      std::function<Eigen::VectorXd(State<ConfigType> const &)> _function);

  std::string name;
  std::string description;
  std::function<Eigen::VectorXd(State<ConfigType> const &)> function;

  Eigen::VectorXd operator()(State<ConfigType> const &state) const;
};

/// StateSampler stores vector valued samples in a matrix
///
/// This just holds a StateSamplingFunction and a Sampler together
template <typename _ConfigType>
class StateSampler {
 public:
  typedef _ConfigType ConfigType;

  /// \brief Sampler constructor
  StateSampler(StateSamplingFunction<ConfigType> _function,
               std::shared_ptr<Sampler> const &_generic_sampler =
                   std::make_shared<Sampler>());

  /// \brief Add a new sample
  void sample(State<ConfigType> const &state);

  /// Get sampling function
  StateSamplingFunction<ConfigType> const &function() const;

  /// Get generic sampler
  std::shared_ptr<Sampler> const &sampler() const;

 private:
  /// Function: [](State<ConfigType> const &state) -> Eigen::VectorXd
  StateSamplingFunction<ConfigType> m_function;

  std::shared_ptr<Sampler> m_sampler;
};

template <typename ConfigType>
std::unique_ptr<StateSampler<ConfigType>> make_conditions_sampler(
    std::string name, std::string description, std::string condition_name);

template <typename ConfigType>
std::unique_ptr<StateSampler<ConfigType>> make_properties_sampler(
    std::string name, std::string description, std::string property_name);

template <typename ConfigType>
std::unique_ptr<StateSampler<ConfigType>> make_configuration_sampler(
    std::string name, std::string description,
    std::function<Eigen::VectorXd(ConfigType const &)> function);

}  // namespace Monte2
}  // namespace CASM

// --- Inline implementations ---

namespace CASM {
namespace Monte2 {

/// \brief StateSamplingFunction constructor
template <typename _ConfigType>
StateSamplingFunction<_ConfigType>::StateSamplingFunction(
    std::string _name, std::string _description,
    std::function<Eigen::VectorXd(State<ConfigType> const &)> _function)
    : name(_name), description(_description), function(_function) {}

/// \brief Take a sample
template <typename _ConfigType>
Eigen::VectorXd StateSamplingFunction<_ConfigType>::operator()(
    State<ConfigType> const &state) const {
  return function(state);
}

/// \brief StateSampler constructor
template <typename _ConfigType>
StateSampler<_ConfigType>::StateSampler(
    StateSamplingFunction<ConfigType> _function,
    std::shared_ptr<Sampler> const &_generic_sampler)
    : m_function(_function), m_sampler(_generic_sampler) {}

/// \brief Add a new sample
template <typename _ConfigType>
void StateSampler<_ConfigType>::sample(State<ConfigType> const &state) {
  m_sampler->push_back(m_function(state));
}

/// Get sampling function
template <typename _ConfigType>
StateSamplingFunction<_ConfigType> const &StateSampler<_ConfigType>::function()
    const {
  return m_function;
}

/// Get generic sampler
template <typename _ConfigType>
std::shared_ptr<Sampler> const &StateSampler<_ConfigType>::sampler() const {
  return m_sampler;
}

template <typename ConfigType>
std::unique_ptr<StateSampler<ConfigType>> make_conditions_sampler(
    std::string name, std::string description, std::string condition_name) {
  auto lambda = [=](State<ConfigType> const &state) {
    return state.conditions.at(condition_name);
  };
  StateSamplingFunction<ConfigType> f(name, description, lambda);
  return notstd::make_unique<StateSampler<ConfigType>>(f);
}

template <typename ConfigType>
std::unique_ptr<StateSampler<ConfigType>> make_properties_sampler(
    std::string name, std::string description, std::string property_name) {
  auto lambda = [=](State<ConfigType> const &state) {
    return state.properties.at(property_name);
  };
  StateSamplingFunction<ConfigType> f(name, description, lambda);
  return notstd::make_unique<StateSampler<ConfigType>>(f);
}

template <typename ConfigType>
std::unique_ptr<StateSampler<ConfigType>> make_configuration_sampler(
    std::string name, std::string description,
    std::function<Eigen::VectorXd(ConfigType const &)> function) {
  auto lambda = [=](State<ConfigType> const &state) {
    return function(state.configuration);
  };
  StateSamplingFunction<ConfigType> f(name, description, lambda);
  return notstd::make_unique<StateSampler<ConfigType>>(f);
}

}  // namespace Monte2
}  // namespace CASM

#endif
