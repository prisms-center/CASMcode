#ifndef CASM_monte2_SamplingParams
#define CASM_monte2_SamplingParams

#include <vector>

#include "casm/monte2/definitions.hh"

namespace CASM {
namespace Monte2 {

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

}  // namespace Monte2
}  // namespace CASM

// --- Inline implementations ---

namespace CASM {
namespace Monte2 {

/// Default constructor
inline SamplingParams::SamplingParams()
    : sample_mode(SAMPLE_MODE::BY_PASS),
      sample_method(SAMPLE_METHOD::LINEAR),
      count_sampling_params({0, 1}),
      time_sampling_params({0., 1.}),
      sampler_names({}),
      sample_trajectory(false) {}

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

}  // namespace Monte2
}  // namespace CASM

#endif
