#include "casm/crystallography/SpeciesAttribute.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/external/Eigen/Core"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
  namespace xtal {

    bool SpeciesAttribute::identical(SpeciesAttribute const &other, double _tol) const {
      return name() == other.name() && almost_equal(value(), other.value(), _tol);
    }
  } // namespace xtal

} // namespace CASM

namespace CASM {
  namespace sym {
    xtal::SpeciesAttribute &apply(const xtal::SymOp &op, xtal::SpeciesAttribute &mutating_attribute) {
      Eigen::MatrixXd symop_matrix_representation =
        mutating_attribute.traits().symop_to_matrix(get_matrix(op), get_translation(op), get_time_reversal(op));
      Eigen::Vector3d new_value = symop_matrix_representation * mutating_attribute.value();
      mutating_attribute.set_value(new_value);
      return mutating_attribute;
    }

    xtal::SpeciesAttribute copy_apply(const xtal::SymOp &op, xtal::SpeciesAttribute attribute) {
      apply(op, attribute);
      return attribute;
    }
  } // namespace sym
} // namespace CASM
