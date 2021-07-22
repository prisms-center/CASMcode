#include "casm/crystallography/SpeciesProperty.hh"

#include "casm/crystallography/SymType.hh"
#include "casm/external/Eigen/Core"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
namespace xtal {

bool SpeciesProperty::identical(SpeciesProperty const &other,
                                double _tol) const {
  return name() == other.name() && almost_equal(value(), other.value(), _tol);
}
}  // namespace xtal

}  // namespace CASM

namespace CASM {
namespace sym {
xtal::SpeciesProperty &apply(const xtal::SymOp &op,
                             xtal::SpeciesProperty &mutating_attribute) {
  Eigen::MatrixXd symop_matrix_representation =
      mutating_attribute.traits().symop_to_matrix(
          get_matrix(op), get_translation(op), get_time_reversal(op));
  mutating_attribute.set_value(symop_matrix_representation *
                               mutating_attribute.value());
  return mutating_attribute;
}

xtal::SpeciesProperty copy_apply(const xtal::SymOp &op,
                                 xtal::SpeciesProperty attribute) {
  apply(op, attribute);
  return attribute;
}
}  // namespace sym
}  // namespace CASM
