#include "casm/crystallography/Site.hh"
#include <casm/crystallography/DoFSet.hh>
#include <casm/crystallography/SymType.hh>
#include <unordered_set>
#include <vector>

namespace CASM {
  namespace xtal {
    bool DoFSetIsEquivalent_f::_traits_match(const DoFSet &other_value) const {
      return m_reference_dofset.traits().name() == other_value.traits().name();
    }

    bool DoFSetIsEquivalent_f::_axis_names_match(const DoFSet &other_value) const {
      return m_reference_dofset.component_names() == other_value.component_names();
    }

    bool DoFSetIsEquivalent_f::_basis_is_equivalent(const DoFSet &other_value) const {
      const auto &reference_basis = m_reference_dofset.basis();
      const auto &other_basis = other_value.basis();

      if(reference_basis.cols() != other_basis.cols())
        return false;

      // If not a square matrix, make sure that column-space _after_basis is identical to that of _before_basis
      if(reference_basis.cols() < other_basis.rows()) {

        // Find rank of augmented matrix. If it is the same as rank of DoF Basis, then
        // the two matrices are similar (and, thus, equivalent)
        Eigen::MatrixXd aug(reference_basis.rows(), 2 * reference_basis.cols());
        aug << reference_basis, other_basis;

        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(aug);
        qr.setThreshold(m_tol);
        if(qr.rank() != m_reference_dofset.dimensions())
          return false;
      }

      return true;
    }

    bool DoFSetIsEquivalent_f::operator()(const DoFSet &other_value) const {
      return this->_traits_match(other_value) && this->_axis_names_match(other_value) && this->_basis_is_equivalent(other_value);
    }

  } // namespace xtal
}

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
