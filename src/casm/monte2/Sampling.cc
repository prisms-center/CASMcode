#include "casm/monte2/Sampling.hh"

namespace CASM {
namespace Monte2 {

/// Lookup table for state sampling functions
std::map<std::string, StateSamplingFunction> &state_sampling_functions() {
  static std::map<std::string, StateSamplingFunction> _state_sampling_functions;
  return _state_sampling_functions;
}

std::unique_ptr<Sampler> make_conditions_sampler(std::string name,
                                                 std::string description,
                                                 std::string condition_name) {
  auto lambda = [=](State const &state) {
    return state.conditions.at(condition_name);
  };
  StateSamplingFunction f(name, description, lambda);
  return notstd::make_unique<Sampler>(f);
}

std::unique_ptr<Sampler> make_properties_sampler(std::string name,
                                                 std::string description,
                                                 std::string property_name) {
  auto lambda = [=](State const &state) {
    return state.properties.at(property_name);
  };
  StateSamplingFunction f(name, description, lambda);
  return notstd::make_unique<Sampler>(f);
}

std::unique_ptr<Sampler> make_configuration_sampler(
    std::string name, std::string description,
    std::function<Eigen::VectorXd(Configuration const &)> function) {
  auto lambda = [=](State const &state) {
    return function(state.configuration);
  };
  StateSamplingFunction f(name, description, lambda);
  return notstd::make_unique<Sampler>(f);
}

/// Default constructor
SamplingParams::SamplingParams()
    : sample_mode(SAMPLE_MODE::BY_PASS),
      sample_method(SAMPLE_METHOD::LINEAR),
      count_sampling_params({0, 1}),
      time_sampling_params({0., 1.}),
      sampler_names({}),
      sample_trajectory(false) {}

}  // namespace Monte2
}  // namespace CASM
