#ifndef CASM_monte2_Sampling
#define CASM_monte2_Sampling

#include <vector>

#include "casm/monte2/State.hh"

namespace CASM {
namespace Monte2 {

struct StateSamplingFunction {
  StateSamplingFunction(
      std::string _name, std::string _description,
      std::function<Eigen::VectorXd(State const &)> _function);

  std::string name;
  std::string description;
  std::function<Eigen::VectorXd(State const &)> function;

  Eigen::VectorXd operator()(State const &state) const;
};

/// Lookup table for state sampling functions
std::map<std::string, StateSamplingFunction> &state_sampling_functions();

/// Sampler stores vector valued samples in a matrix
class Sampler {
 public:
  /// \brief Sampler constructor
  Sampler(StateSamplingFunction _function,
          CountType _capacity_increment = 1000);

  /// \brief Add a new sample
  void sample(State const &state);

  /// \brief Add a new sample
  void push_back(Eigen::VectorXd const &vector);

  /// \brief Set all values directly
  void set_values(Eigen::MatrixXd const &values);

  /// Clear values - preserves n_components, set n_samples to 0
  void clear();

  /// Non-conservative resize - sets n_samples to 0
  void resize(Index n_components);

  /// Conservative resize, to increase capacity for more samples
  void set_sample_capacity(CountType sample_capacity);

  /// Set capacity increment (used when push_back requires more capacity)
  void set_capacity_increment(CountType _capacity_increment);

  /// Get sampler name
  std::string name() const;

  /// Get sampling function
  StateSamplingFunction const &function() const;

  /// Number of components (vector size) of samples
  Index n_components() const;

  /// Current number of samples taken
  CountType n_samples() const;

  /// Current sample capacity
  CountType sample_capacity() const;

  /// \brief Get sampled values as a matrix
  ///
  /// Stores individual vector samples in rows,
  /// use columns to check convergence of individual components
  Eigen::Block<const Eigen::MatrixXd> values() const;

  /// \brief Get all samples of a particular component (a column of `values()`)
  Eigen::Block<const Eigen::MatrixXd> component(Index component_index) const;

  /// \brief Get a sample (a row of `values()`)
  Eigen::Block<const Eigen::MatrixXd> sample(CountType sample_index) const;

 private:
  /// Function: [](State const &state) -> Eigen::VectorXd
  StateSamplingFunction m_function;

  /// Current number of samples taken
  Index m_n_samples;

  /// Used when push_back requires more capacity
  CountType m_capacity_increment;

  /// Stores individual samples in rows.
  /// Total number of rows is `sample_capacity` which to avoid constant re-
  /// sizing may be greater than the current number of samples taken,
  /// `m_n_samples`
  /// Use columns to check convergence of individual components.
  Eigen::MatrixXd m_values;
};

struct SamplerComponent {
  /// Name (i.e. "comp_n(Mg)", "corr(0)", etc.)
  std::string name = "";

  /// Sampler name (i.e. "comp_n", "corr", etc.)
  std::string sampler_name = "";

  /// Sampler component index (i.e. 0, 1, etc.)
  Index component_index = 0;

  bool operator<(SamplerComponent const &other) const;
};

std::map<std::string, Sampler>::const_iterator find_and_validate(
    SamplerComponent const &key,
    std::map<std::string, Sampler> const &samplers);

std::unique_ptr<Sampler> make_conditions_sampler(std::string name,
                                                 std::string description,
                                                 std::string condition_name);

std::unique_ptr<Sampler> make_properties_sampler(std::string name,
                                                 std::string description,
                                                 std::string property_name);

std::unique_ptr<Sampler> make_configuration_sampler(
    std::string name, std::string description,
    std::function<Eigen::VectorXd(Configuration const &)> function);

/// How often to sample runs
enum class SAMPLE_MODE { BY_STEP, BY_PASS, BY_TIME };

/// How to sample by time
enum class SAMPLE_METHOD { LINEAR, LOG };

/// What to sample and how
struct SamplingParams {
  /// Default constructor
  SamplingParams();

  /// \brief Sample by step, pass, or time
  ///
  /// Default=SAMPLE_MODE::BY_PASS
  SAMPLE_MODE sample_mode;

  /// \brief Sample linearly or logarithmically
  ///
  /// Default=SAMPLE_METHOD::LINEAR
  SAMPLE_METHOD sample_method;

  // --- Used for sampling BY_STEP or BY_PASS ---

  /// \brief Parameters for determining when samples are taken by count
  ///
  /// Linear sampling: count = a + b*n
  /// Log sampling: count = a + b^(n-c)
  ///   where n=sample index
  ///
  /// Default={0, 1} (sample period = 1)
  std::vector<CountType> count_sampling_params;

  // --- Used for sampling BY_TIME ---

  /// \brief Parameters for determining when samples are taken by time
  ///
  /// Linear sampling: time = a + b*n
  /// Log sampling: time = a + b^(n-c)
  ///   where n=sample index, starting from 0
  ///
  /// Default={0., 1.} (sample period = 1.)
  std::vector<TimeType> time_sampling_params;

  /// \brief What to sample
  ///
  /// Default={}
  std::vector<std::string> sampler_names;

  /// \brief Whether to sample the complete configuration
  ///
  /// Default=false
  bool sample_trajectory;
};

/// The count when a particular sample is due
CountType sample_count(SAMPLE_METHOD sample_method,
                       std::vector<CountType> const &count_sampling_params,
                       CountType sample_index);

/// The time when a particular sample is due
TimeType sample_time(SAMPLE_METHOD sample_method,
                     std::vector<TimeType> const &time_sampling_params,
                     CountType sample_index);

// /// Tracks step, pass, and / or time
// class Counter {
//  public:
//   Counter(CountType _steps_per_pass)
//       : pass(0), step(0), steps_per_pass(_steps_per_pass), time(0.0) {}
//
//   /// One pass is equal to the number of sites with variable degrees of
//   freedom
//   /// in the supercell
//   CountType pass;
//
//   /// Number of steps taken during the current pass; gets reset to zero when
//   it
//   /// reaches the number of steps per pass
//   CountType step;
//
//   /// One pass is equal to the number of sites with variable degrees of
//   freedom
//   /// in the supercell
//   CountType steps_per_pass;
//
//   /// Simulated time
//   TimeType time;
//
//   /// Increment step, then update pass and step based on steps_per_pass
//   Counter &increment_step();
//
//   /// Increment time
//   Counter &increment_time(TimeType time_increment);
// };

// Sampled data from a single run (constant conditions)
struct SampledData {
  /// Map of <sampler name>:<sampler>
  std::map<std::string, Sampler> samplers;

  /// Samples of the complete configuration
  std::vector<Configuration> trajectory;

  /// Vector of counts (could be pass or step) when a sample occurred
  std::vector<CountType> count;

  /// Vector of times when a sample occurred
  std::vector<TimeType> time;
};

/// \brief Last sampled_data.count value (else 0)
CountType get_count(SampledData const &sampled_data);

/// \brief Last sampled_data.time value (else 0)
TimeType get_time(SampledData const &sampled_data);

/// \brief Get Sampler::n_samples() value (assumes same for all) (else 0)
CountType get_n_samples(std::map<std::string, Sampler> const &samplers);

/// \brief Get Sampler::n_samples() value (assumes same for all) (else 0)
CountType get_n_samples(SampledData const &sampled_data);

}  // namespace Monte2
}  // namespace CASM

// --- Inline implementations ---

namespace CASM {
namespace Monte2 {

inline StateSamplingFunction::StateSamplingFunction(
    std::string _name, std::string _description,
    std::function<Eigen::VectorXd(State const &)> _function)
    : name(_name), description(_description), function(_function) {}

inline Eigen::VectorXd StateSamplingFunction::operator()(
    State const &state) const {
  return function(state);
}

/// \brief Sampler constructor
///
/// \param _name Name of value that will be sampled
/// \param _description Extended description of sampler, to be printed with
///     help output
/// \param _capacity_increment How much to resize the underlying matrix by
///     whenever space runs out.
inline Sampler::Sampler(StateSamplingFunction _function,
                        CountType _capacity_increment)
    : m_function(_function),
      m_n_samples(0),
      m_capacity_increment(_capacity_increment) {}

/// \brief Add a new sample
inline void Sampler::sample(State const &state) {
  this->push_back(m_function(state));
}

/// \brief Add a new sample
inline void Sampler::push_back(Eigen::VectorXd const &vector) {
  if (n_samples() == sample_capacity()) {
    set_sample_capacity(sample_capacity() + m_capacity_increment);
  }
  m_values.row(m_n_samples) = vector;
  ++m_n_samples;
}

/// \brief Set all values directly
inline void Sampler::set_values(Eigen::MatrixXd const &values) {
  m_values = values;
  m_n_samples = m_values.rows();
}

/// Clear values - preserves n_components, set n_samples to 0
inline void Sampler::clear() { resize(n_components()); }

/// Non-conservative resize - sets n_samples to 0
inline void Sampler::resize(Index n_components) {
  m_values.resize(m_capacity_increment, n_components);
  m_n_samples = 0;
}

/// Conservative resize, to increase capacity for more samples
inline void Sampler::set_sample_capacity(CountType sample_capacity) {
  m_values.conservativeResize(sample_capacity, Eigen::NoChange_t());
}

/// Set capacity increment (used when push_back requires more capacity)
inline void Sampler::set_capacity_increment(CountType _capacity_increment) {
  m_capacity_increment = _capacity_increment;
}

/// Get sampler name
inline std::string Sampler::name() const { return m_function.name; }

/// Get sampling function
inline StateSamplingFunction const &Sampler::function() const {
  return m_function;
}

/// Number of components (vector size) of samples
inline Index Sampler::n_components() const { return m_values.cols(); }

/// Current number of samples taken
inline CountType Sampler::n_samples() const { return m_n_samples; }

/// Current sample capacity
inline CountType Sampler::sample_capacity() const { return m_values.rows(); }

/// \brief Get sampled values as a matrix
///
/// Stores individual vector samples in rows,
/// use columns to check convergence of individual components
inline Eigen::Block<const Eigen::MatrixXd> Sampler::values() const {
  return m_values.block(0, 0, m_n_samples, m_values.cols());
}

/// \brief Get all samples of a particular component (a column of `values()`)
inline Eigen::Block<const Eigen::MatrixXd> Sampler::component(
    Index component_index) const {
  return m_values.block(0, component_index, m_n_samples, 1);
}

/// \brief Get a sample (a row of `values()`)
inline Eigen::Block<const Eigen::MatrixXd> Sampler::sample(
    CountType sample_index) const {
  return m_values.block(sample_index, 0, 1, m_values.cols());
}

inline bool SamplerComponent::operator<(SamplerComponent const &other) const {
  if (this->sampler_name == other.sampler_name) {
    return this->component_index < other.component_index;
  }
  return this->component_index < other.component_index;
}

inline std::map<std::string, Sampler>::const_iterator find_and_validate(
    SamplerComponent const &key,
    std::map<std::string, Sampler> const &samplers) {
  // find and validate sampler name && component index
  auto sampler_it = samplers.find(key.sampler_name);
  if (sampler_it == samplers.end()) {
    std::stringstream msg;
    msg << "Error in equilibration_check: Sampler '" << key.sampler_name
        << "' not found." << std::endl;
    throw std::runtime_error(msg.str());
  }
  if (key.component_index >= sampler_it->second.n_components()) {
    std::stringstream msg;
    msg << "Error in equilibration_check: Requested component index "
        << key.component_index << ", but '" << key.sampler_name << "' has "
        << sampler_it->second.n_components() << "components." << std::endl;
    throw std::runtime_error(msg.str());
  }
  return sampler_it;
}

inline CountType sample_count(
    SAMPLE_METHOD sample_method,
    std::vector<CountType> const &count_sampling_params,
    CountType sample_index) {
  if (sample_method == SAMPLE_METHOD::LINEAR) {
    auto &a = count_sampling_params[0];
    auto &b = count_sampling_params[1];
    auto &n = sample_index;
    return a + b * n;
  } else /* sample_method == SAMPLE_METHOD::LOG */ {
    auto &a = count_sampling_params[0];
    auto &b = count_sampling_params[1];
    auto &c = count_sampling_params[2];
    auto &n = sample_index;
    return a + pow(b, n - c);
  }
}

inline TimeType sample_time(SAMPLE_METHOD sample_method,
                            std::vector<TimeType> const &time_sampling_params,
                            CountType sample_index) {
  if (sample_method == SAMPLE_METHOD::LINEAR) {
    auto &a = time_sampling_params[0];
    auto &b = time_sampling_params[1];
    auto &n = sample_index;
    return a + b * n;
  } else /* sample_method == SAMPLE_METHOD::LOG */ {
    auto &a = time_sampling_params[0];
    auto &b = time_sampling_params[1];
    auto &c = time_sampling_params[2];
    auto &n = sample_index;
    return a + pow(b, n - c);
  }
}

// inline Counter &Counter::increment_step() {
//   ++step;
//   if (step == steps_per_pass) {
//     ++pass;
//     step = 0;
//   }
//
//   return *this;
// }
//
// inline Counter &Counter::increment_time(TimeType time_increment) {
//   time += time_increment;
//   return *this;
// }

inline CountType get_count(SampledData const &sampled_data) {
  if (sampled_data.count.size()) {
    return sampled_data.count.back();
  }
  return CountType(0);
}

inline TimeType get_time(SampledData const &sampled_data) {
  if (sampled_data.time.size()) {
    return sampled_data.time.back();
  }
  return TimeType(0.0);
}

inline CountType get_n_samples(std::map<std::string, Sampler> const &samplers) {
  if (samplers.size()) {
    return samplers.begin()->second.n_samples();
  }
  return CountType(0);
}

inline CountType get_n_samples(SampledData const &sampled_data) {
  return get_n_samples(sampled_data.samplers);
}

}  // namespace Monte2
}  // namespace CASM

#endif
