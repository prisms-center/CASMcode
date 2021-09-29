#ifndef CASM_monte2_CutoffCheck
#define CASM_monte2_CutoffCheck

#include <optional>

#include "casm/monte2/definitions.hh"

namespace CASM {
namespace Monte2 {

// --- Cutoff checking ---

/// \brief Completion check parameters that don't depend on the sampled values
struct CutoffCheckParams {
  // --- A calculation does not stop before all minimums are met ---

  std::optional<CountType> min_count;
  std::optional<TimeType> min_time;
  std::optional<CountType> min_sample;

  // --- A calculation does stop when any maximum is met ---

  std::optional<CountType> max_count;
  std::optional<TimeType> max_time;
  std::optional<CountType> max_sample;
};

bool all_minimums_met(CutoffCheckParams const &cutoff_params,
                      std::optional<CountType> count,
                      std::optional<TimeType> time, CountType n_samples);

bool any_maximum_met(CutoffCheckParams const &cutoff_params,
                     std::optional<CountType> count,
                     std::optional<TimeType> time, CountType n_samples);

// --- inline definitions ---

inline bool all_minimums_met(CutoffCheckParams const &cutoff_params,
                             std::optional<CountType> count,
                             std::optional<TimeType> time,
                             CountType n_samples) {
  auto const &p = cutoff_params;

  if (p.min_sample.has_value() && n_samples < p.min_sample.value()) {
    return false;
  }

  if (p.min_count.has_value() && count.has_value() &&
      count.value() < p.min_count.value()) {
    return false;
  }

  if (p.min_time.has_value() && time.has_value() &&
      time.value() < p.min_time.value()) {
    return false;
  }

  return true;
}

inline bool any_maximum_met(CutoffCheckParams const &cutoff_params,
                            std::optional<CountType> count,
                            std::optional<TimeType> time, CountType n_samples) {
  auto const &p = cutoff_params;

  if (p.max_sample.has_value() && n_samples >= p.max_sample.value()) {
    return true;
  }

  if (p.max_count.has_value() && count.has_value() &&
      count.value() >= p.max_count.value()) {
    return true;
  }

  if (p.max_time.has_value() && time.has_value() &&
      time.value() >= p.max_time.value()) {
    return true;
  }

  return false;
}

}  // namespace Monte2
}  // namespace CASM

#endif
