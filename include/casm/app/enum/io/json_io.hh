#ifndef CASM_app_enum_json_io
#define CASM_app_enum_json_io

#include <string>
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {

  struct KwargsParser;

  /// Parse a number or array value and return a vector
  Eigen::VectorXd parse_vector_from_number_or_array(
    KwargsParser &parser,
    std::string attribute_name,
    Index dimension,
    Eigen::VectorXd const *default_value = nullptr);

  /// Parse "axes", "min", "max", and "increment" for continuous DoF enumeration
  void parse_dof_space_axes(KwargsParser &parser,
                            Eigen::MatrixXd &axes,
                            Eigen::VectorXd &min_val,
                            Eigen::VectorXd &max_val,
                            Eigen::VectorXd &inc_val,
                            Index dof_space_dimension);
}

#endif
