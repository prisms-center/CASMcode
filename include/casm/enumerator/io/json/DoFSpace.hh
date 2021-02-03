#ifndef CASM_enumerator_io_json_DoFSpace
#define CASM_enumerator_io_json_DoFSpace

#include <string>

#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"
#include "casm/symmetry/SymRepTools.hh"

namespace CASM {

class DoFSpace;
template <typename T>
class InputParser;
class Structure;
template <typename T>
struct jsonConstructor;
template <typename T>
struct jsonMake;
class jsonParser;

// jsonParser &to_json(DoFSpace const &dofspace, jsonParser &json);

void from_json(DoFSpace &dofspace, jsonParser const &json,
               std::shared_ptr<Structure const> const &shared_prim);

jsonParser &to_json(
    DoFSpace const &dofspace, jsonParser &json,
    std::optional<std::string> const &identifier = std::nullopt,
    std::optional<ConfigEnumInput> const &input_state = std::nullopt,
    std::optional<VectorSpaceSymReport> const &sym_report = std::nullopt);

template <>
struct jsonConstructor<DoFSpace> {
  static DoFSpace from_json(
      jsonParser const &json,
      std::shared_ptr<Structure const> const &shared_prim);
};

template <>
struct jsonMake<DoFSpace> {
  static std::unique_ptr<DoFSpace> make_from_json(
      jsonParser const &json,
      std::shared_ptr<Structure const> const &shared_prim);
};

/// Data structure used for continuous DoF enumeration IO
struct AxesCounterParams {
  // DoFSpace axes (column vector matrix, may be a subspace)

  Eigen::MatrixXd axes;

  // Counter min, max, and increment values (size must match axes.cols())

  Eigen::VectorXd min_val;
  Eigen::VectorXd max_val;
  Eigen::VectorXd inc_val;
};

/// Parse DoFSpace subspace from  "axes" and normal coordinate grid counter from
/// "min", "max", and "increment"
void parse(InputParser<AxesCounterParams> &parser, Index dof_space_dimension);

/// Adds homogeneous mode info to the json file
void add_homogeneous_mode_info(jsonParser &json, DoFSpace const &dofspace,
                               ConfigEnumInput const &input_state);
}  // namespace CASM

#endif
