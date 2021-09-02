#ifndef CASM_clexulator_ConfigDoFValues
#define CASM_clexulator_ConfigDoFValues
#include <map>
#include <string>

#include "casm/global/eigen.hh"

namespace CASM {
namespace clexulator {

/// \brief Data structure holding configuration DoF values
///
/// Notes:
/// - Site DoF values are organized in blocks by sublattice, and ordered within
///   a sublattice block by unit cell index
/// - For use by Clexulator, DoF values should be expressed in terms of the
///   `prim basis`.
/// - For input / output, conversion to structure, etc. the DoF values should be
///   expressed in terms of the `standard basis`.
/// - See `ConfigDoF` for full documentation on layout.
struct ConfigDoFValues {
  /// \brief Occupation values (shape=(N_sites,1))
  Eigen::VectorXi occupation;

  /// \brief Local continuous DoF values
  ///
  /// For use by clexulator, expected shape is (max(DoFSet::dim()), N_sites)
  /// - The number of rows is the maximum of the individual sublattice DoFSet
  /// dimensions.
  ///
  /// For the standard basis, the expected shape is (standard_dim, N_sites)
  std::map<std::string, Eigen::MatrixXd> local_dof_values;

  /// \brief Global continuous DoF values
  ///
  /// For use by clexulator, expected shape is (DoFSet::dim(), 1)
  ///
  /// For the standard basis, the expected shape is (standard_dim, 1)
  std::map<std::string, Eigen::VectorXd> global_dof_values;
};

}  // namespace clexulator
}  // namespace CASM

#endif
