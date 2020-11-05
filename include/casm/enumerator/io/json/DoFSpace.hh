#ifndef CASM_enumerator_io_json_DoFSpace
#define CASM_enumerator_io_json_DoFSpace

#include <string>
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"


namespace CASM {

  struct DoFSpace;
  template<typename T> class InputParser;

  jsonParser &to_json(DoFSpace const &dofspace, jsonParser &json, std::string name);

  /// Data structure used for continuous DoF enumeration IO
  struct AxesCounterParams {

    // DoFSpace axes (column vector matrix, may be a subspace)

    Eigen::MatrixXd axes;

    // Counter min, max, and increment values (size must match axes.cols())

    Eigen::VectorXd min_val;
    Eigen::VectorXd max_val;
    Eigen::VectorXd inc_val;
  };

  /// Parse DoFSpace subspace from  "axes" and normal coordinate grid counter from "min", "max", and "increment"
  void parse(InputParser<AxesCounterParams> &parser,
             Index dof_space_dimension);
}

#endif
