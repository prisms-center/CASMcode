#include <casm/crystallography/DoFSet.hh>
#include <casm/crystallography/SymType.hh>

namespace CASM {
  namespace sym {
    /// \brief Copy and apply SymOp to a DoFSet
    xtal::DoFSet copy_apply(const xtal::SymOp &op, const xtal::DoFSet &dof) {
      Eigen::Matrix3d transformation_matrix = dof.traits().symop_to_matrix(get_matrix(op), get_translation(op), get_time_reversal(op));
      Eigen::Matrix3d new_basis = transformation_matrix * dof.basis();
      return xtal::DoFSet(dof.traits(), dof.component_names(), new_basis);
    }
  } // namespace sym
} // namespace CASM
