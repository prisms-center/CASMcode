#ifndef CASM_monte2_Results
#define CASM_monte2_Results

#include "casm/clex/Configuration.hh"
#include "casm/monte2/CompletionCheck.hh"
#include "casm/monte2/Definitions.hh"
#include "casm/monte2/Sampling.hh"

namespace CASM {
namespace Monte2 {

/// name -> (mean, precision, etc. as calculated by `summarize` method)
typedef std::map<std::string, VectorValueMap> VectorStatisticsMap;

/// \brief Holds detailed results of a calculation at a single set of
/// conditions, including all of the sampled values
struct Results {
 public:
  /// Initial configuration
  std::optional<Configuration> initial_configuration;

  /// Final configuration
  std::optional<Configuration> final_configuration;

  /// Conditions
  VectorValueMap conditions;

  /// All sampled data
  std::shared_ptr<SampledData> sampled_data;

  /// Completion check results
  std::optional<CompletionCheckResults> completion_check_results;

  /// Statistics from sampled values, calculated from samples in the range
  /// [begin, end), where begin=N_samples_for_all_to_equilibrate, and
  /// end=(N_samples_for_all_to_equilibrate + N_samples_for_statistics)
  ///
  /// Statistics calculated are determined by the Monte Carlo method, but
  /// typically include "mean", "precision" among other values. The method
  /// description should define the values included here.
  std::map<std::string, std::any> statistics;
};

}  // namespace Monte2
}  // namespace CASM

#endif
