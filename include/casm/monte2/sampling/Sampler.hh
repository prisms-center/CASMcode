#ifndef CASM_monte2_Sampler
#define CASM_monte2_Sampler

#include <vector>

#include "casm/global/eigen.hh"
#include "casm/monte2/definitions.hh"

namespace CASM {
namespace Monte2 {

/// \brief Sampler stores vector valued samples in a matrix
class Sampler {
 public:
  /// \brief Sampler constructor
  Sampler(CountType _capacity_increment = 1000);

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

std::map<std::string, std::shared_ptr<Sampler>>::const_iterator
find_and_validate(
    SamplerComponent const &key,
    std::map<std::string, std::shared_ptr<Sampler>> const &samplers);

/// \brief Get Sampler::n_samples() value (assumes same for all)
/// (else 0)
CountType get_n_samples(
    std::map<std::string, std::shared_ptr<Sampler>> const &samplers);

}  // namespace Monte2
}  // namespace CASM

// --- Inline implementations ---

namespace CASM {
namespace Monte2 {

/// \brief Sampler constructor
///
/// \param _name Name of value that will be sampled
/// \param _description Extended description of sampler, to be printed with
///     help output
/// \param _capacity_increment How much to resize the underlying matrix by
///     whenever space runs out.
inline Sampler::Sampler(CountType _capacity_increment)
    : m_n_samples(0), m_capacity_increment(_capacity_increment) {}

/// \brief Add a new sample
inline void Sampler::push_back(Eigen::VectorXd const &vector) {
  if (n_samples() == 0) {
    resize(vector.size());
  }
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

inline std::map<std::string, std::shared_ptr<Sampler>>::const_iterator
find_and_validate(
    SamplerComponent const &key,
    std::map<std::string, std::shared_ptr<Sampler>> const &samplers) {
  // find and validate sampler name && component index
  auto sampler_it = samplers.find(key.sampler_name);
  if (sampler_it == samplers.end()) {
    std::stringstream msg;
    msg << "Error in equilibration_check: Sampler '" << key.sampler_name
        << "' not found." << std::endl;
    throw std::runtime_error(msg.str());
  }
  if (key.component_index >= sampler_it->second->n_components()) {
    std::stringstream msg;
    msg << "Error in find_and_validate: Requested component index "
        << key.component_index << ", but '" << key.sampler_name << "' has "
        << sampler_it->second->n_components() << "components." << std::endl;
    throw std::runtime_error(msg.str());
  }
  return sampler_it;
}

inline CountType get_n_samples(
    std::map<std::string, std::shared_ptr<Sampler>> const &samplers) {
  if (samplers.size()) {
    return samplers.begin()->second->n_samples();
  }
  return CountType(0);
}

}  // namespace Monte2
}  // namespace CASM

#endif
