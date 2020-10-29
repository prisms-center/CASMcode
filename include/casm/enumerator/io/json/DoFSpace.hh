#ifndef CASM_enumerator_io_json_DoFSpace
#define CASM_enumerator_io_json_DoFSpace

#include <string>
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"


namespace CASM {

  struct DoFSpace;
  struct KwargsParser;

  jsonParser &to_json(DoFSpace const &dofspace, jsonParser &json, std::string name);

  /// Parse DoF space "axes"
  void parse_dof_space_axes(KwargsParser &parser,
                            Eigen::MatrixXd &axes,
                            Index dof_space_dimension);

  /// Parse DoF space "axes" along with "min", "max", and "increment"
  void parse_dof_space_axes(KwargsParser &parser,
                            Eigen::MatrixXd &axes,
                            Eigen::VectorXd &min_val,
                            Eigen::VectorXd &max_val,
                            Eigen::VectorXd &inc_val,
                            Index dof_space_dimension);
}

#endif
