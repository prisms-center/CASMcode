#ifndef CASM_monte2_SampledData
#define CASM_monte2_SampledData

#include <vector>

#include "casm/global/eigen.hh"
#include "casm/monte2/sampling/Sampler.hh"

namespace CASM {
namespace Monte2 {

// Sampled data from a single run (constant conditions)
struct SampledData {
  /// Map of <sampler name>:<sampler>
  std::map<std::string, std::shared_ptr<Sampler>> samplers;

  /// Vector of counts (could be pass or step) when a sample occurred
  std::vector<CountType> count;

  /// Vector of times when a sample occurred
  std::vector<TimeType> time;
};

/// \brief Last sampled_data.count value (else 0)
std::optional<CountType> get_count(SampledData const &sampled_data);

/// \brief Last sampled_data.time value (else 0)
std::optional<TimeType> get_time(SampledData const &sampled_data);

/// \brief Get Sampler::n_samples() value (assumes same for all)
/// (else 0)
CountType get_n_samples(SampledData const &sampled_data);

}  // namespace Monte2
}  // namespace CASM

// --- Inline implementations ---

namespace CASM {
namespace Monte2 {

inline std::optional<CountType> get_count(SampledData const &sampled_data) {
  if (sampled_data.count.size()) {
    return sampled_data.count.back();
  }
  return std::nullopt;
}

inline std::optional<TimeType> get_time(SampledData const &sampled_data) {
  if (sampled_data.time.size()) {
    return sampled_data.time.back();
  }
  return std::nullopt;
}

inline CountType get_n_samples(SampledData const &sampled_data) {
  return get_n_samples(sampled_data.samplers);
}

}  // namespace Monte2
}  // namespace CASM

#endif
