#ifndef CASM_monte2_StateGenerator
#define CASM_monte2_StateGenerator

#include <vector>

#include "casm/monte2/definitions.hh"
#include "casm/monte2/state/State.hh"

namespace CASM {
namespace Monte2 {

/// \brief A StateGenerator generates initial states for a series of Monte
///     Carlo calculations
///
/// Notes:
/// - The template parameter _RunInfoType is specified by a particular Monte
///   Carlo method implementation.
/// - _RunInfoType allows customization of what
///   information is provided to a particular state generation method. In the
///   basic case, it will be the final state for each run. Templating allows
///   support for more complex cases where the next state could be generated
///   based on the sampled data generated during previous runs.
template <typename _ConfigType, typename _RunInfoType>
class StateGenerator {
 public:
  typedef _ConfigType ConfigType;
  typedef _RunInfoType RunInfoType;

  /// \brief Check if calculations are complete, using info from all finished
  ///     runs
  virtual bool is_complete(std::vector<RunInfoType> const &run_info) = 0;

  /// \brief Generate the next initial state, using info from all finished
  ///     runs
  virtual State<ConfigType> next_state(
      std::vector<RunInfoType> const &run_info) = 0;
};

}  // namespace Monte2
}  // namespace CASM

#endif
