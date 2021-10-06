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
  /// \brief Sampler constructor - default component names
  Sampler(Index n_components, CountType _capacity_increment = 1000);

  /// \brief Sampler constructor - custom component names
  Sampler(std::vector<std::string> const &_component_names,
          CountType _capacity_increment = 1000);

  /// \brief Add a new sample
  void push_back(Eigen::VectorXd const &vector);

  /// \brief Set all values directly
  void set_values(Eigen::MatrixXd const &values);

  /// Clear values - preserves n_components, set n_samples to 0
  void clear();

  /// Conservative resize, to increase capacity for more samples
  void set_sample_capacity(CountType sample_capacity);

  /// Set capacity increment (used when push_back requires more capacity)
  void set_capacity_increment(CountType _capacity_increment);

  std::vector<std::string> component_names()

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
  /// Size of vectors to be sampled
  Index m_n_components;

  /// Names to use for components. Size == m_n_components
  std::vector<std::string> m_component_names;

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
  /// \brief Constructor
  ///
  /// \param _sampler_name Name of Sampler
  /// \param _component_index Index into sampler output vector
  SamplerComponent(std::string _sampler_name, Index _component_index)
      : sampler_name(_sampler_name), component_index(_component_index) {}

  /// Sampler name (i.e. "comp_n", "corr", etc.)
  std::string sampler_name = "";

  /// Sampler component index (i.e. 0, 1, etc.)
  Index component_index = 0;

  bool operator<(SamplerComponent const &other) const;
};

struct SamplerConvergenceParams {
  SamplerConvergenceParams(double _precision, double _confidence = 0.95)
      : precision(_precision), confidence(_confidence) {}

  double precision;

  double confidence;
};

/// \brief Helper for compact construction of sampler convergence params
///
/// Usage:
/// - This class is intended as a temporary intermediate constructed by the
/// `converge` function. See `converge` documentation for intended usage.
struct SamplerConvergenceParamsConstructor {
  /// \brief Constructor
  ///
  /// \param _sampler_name Sampler name
  /// \param _sampler Sampler reference
  ///
  /// Note:
  /// - Constructs `values` to include convergence parameters for all
  ///   components of the specified sampler, with initial values precision =
  ///   `std::numeric_limits<double>::infinity()`, confidence = `0.95`.
  SamplerConvergenceParamsConstructor(std::string _sampler_name,
                                      Sampler const &_sampler)
      : sampler_name(_sampler_name), sampler(_sampler) {
    for (Index i = 0; i < sampler.n_components(); ++i) {
      values.emplace(
          SamplerComponent(sampler_name, i),
          SamplerConvergenceParams(std::numeric_limits<double>::infinity()));
    }
  }

  /// \brief Select only the specified component - by index
  SamplerConvergenceParamsConstructor &component(Index component_index) {
    if (component_index >= sampler.component_names().size()) {
      std::stringstream msg;
      msg << "Error constructing sampler convergence parameters: Component "
             "index '"
          << component_name << "' out of range for sampler '" << sampler_name
          << "'";
      throw std::runtime_error(msg.str());
    }
    SamplerComponent component(sampler_name, component_index);
    SamplerConvergenceParams chosen = values.at(component);
    values.clear();
    values.emplace(component, chosen);
    return *this;
  }

  /// \brief Select only the specified component - by name
  SamplerConvergenceParamsConstructor &component(Index component_name) {
    auto begin = sampler.component_names().begin();
    auto end = sampler.component_names().end();
    auto it = std::find(begin, end, component_name);
    if (it == end) {
      std::stringstream msg;
      msg << "Error constructing sampler convergence parameters: Cannot find "
             "component '"
          << component_name << "' for sampler '" << sampler_name << "'";
      throw std::runtime_error(msg.str());
    }
    Index component_index = std::distance(begin, it);
    SamplerComponent component(sampler_name, component_index);
    SamplerConvergenceParams chosen = values.at(component);
    values.clear();
    values.emplace(component, chosen);
    return *this;
  }

  /// \brief Set the requested convergence precision for all selected components
  SamplerConvergenceParamsConstructor &precision(double _precision) {
    for (auto &value : values) {
      value.second.precision = _precision;
    }
    return *this;
  }

  /// \brief Set the requested convergence confidence level
  SamplerConvergenceParamsConstructor &confidence(double _confidence) {
    for (auto &value : values) {
      value.second.confidence = _confidence;
    }
    return *this;
  }

  std::string sampler_name;
  Sampler const &sampler;
  std::map<SamplerComponent, SamplerConvergenceParams> values;
};

std::vector<std::map<SamplerComponent, SamplerConvergence>> converge

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

/// \brief Sampler constructor - default component names
///
/// \param _n_components Size of vectors to be sampled
/// \param _capacity_increment How much to resize the underlying matrix by
///     whenever space runs out.
///
/// Notes:
/// - Components are given default names ["0", "1", "2", ...)
inline Sampler::Sampler(Index _n_components, CountType _capacity_increment)
    : m_n_components(_n_components),
      m_n_samples(0),
      m_capacity_increment(_capacity_increment) {
  for (Index i = 0; i < m_n_components; ++i) {
    m_component_names.push_back(std::to_string(i));
  }
}

/// \brief Sampler constructor - custom component names
///
/// \param _component_names Names to give to each sampled vector element
/// \param _capacity_increment How much to resize the underlying matrix by
///     whenever space runs out.
inline Sampler::Sampler(std::vector<std::string> const &_component_names,
                        CountType _capacity_increment)
    : m_n_components(_component_names.size()),
      m_component_names(_component_names),
      m_n_samples(0),
      m_capacity_increment(_capacity_increment) {}

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
inline void Sampler::clear() {
  m_values.resize(m_capacity_increment, m_n_components);
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
inline Index Sampler::n_components() const { return m_n_components; }

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
  return this->sampler_name < other.sampler_name;
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
