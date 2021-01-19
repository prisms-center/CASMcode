#include "casm/crystallography/IntegralCoordinateWithin.hh"

#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/external/Eigen/Core"
#include "casm/global/eigen.hh"

namespace CASM {
namespace xtal {

void IntegralCoordinateWithin_f::_throw_if_bad_transformation_matrix(
    const matrix_type &transformation_matrix) {
  if (transformation_matrix.determinant() == 0) {
    throw std::runtime_error(
        "The transformation matrix that converts the tiling unit to the "
        "superlattice is singular, and therefore not valid.");
  }
  return;
}

IntegralCoordinateWithin_f::vector_type IntegralCoordinateWithin_f::operator()(
    const vector_type &ijk) const {
  vector_type vec2 = this->m_transformation_matrix_adjugate * ijk;
  long vol = this->m_transformation_matrix
                 .determinant();  // this volume is signed, could be negative,
                                  // that's fine.
  vec2[0] = ((vec2[0] % vol) + std::abs(vol)) % vol;
  vec2[1] = ((vec2[1] % vol) + std::abs(vol)) % vol;
  vec2[2] = ((vec2[2] % vol) + std::abs(vol)) % vol;
  return (this->m_transformation_matrix * vec2) / vol;
}

UnitCellCoord IntegralCoordinateWithin_f::operator()(
    const UnitCellCoord &bijk) const {
  return UnitCellCoord(bijk.sublattice(), this->operator()(bijk.unitcell()));
}

//********************************************************************************************************************************//

namespace impl {
OrderedLatticePointGenerator::OrderedLatticePointGenerator(
    const matrix_type &transformation_matrix)
    : m_bring_within_f(transformation_matrix),
      m_total_lattice_points(std::abs(transformation_matrix.determinant())) {
  smith_normal_form(transformation_matrix, this->m_smith_normal_U,
                    this->m_smith_normal_S, this->m_smith_normal_V);

  vector_type S_diagonal = this->m_smith_normal_S.diagonal();

  m_stride[0] = S_diagonal[0];
  m_stride[1] = S_diagonal[0] * S_diagonal[1];

  assert(transformation_matrix ==
         m_smith_normal_U * m_smith_normal_S * m_smith_normal_V);
}

OrderedLatticePointGenerator::vector_type
OrderedLatticePointGenerator::_normalize_lattice_point(
    const vector_type &mnp) const {
  vector_type ijk = m_smith_normal_U * mnp;
  return m_bring_within_f(ijk);
}

OrderedLatticePointGenerator::vector_type
OrderedLatticePointGenerator::_make_smith_normal_form_lattice_point(
    Index ix) const {
  vector_type mnp((ix % m_stride[1]) % m_stride[0],
                  (ix % m_stride[1]) / m_stride[0], ix / m_stride[1]);
  return mnp;
}

OrderedLatticePointGenerator::vector_type
OrderedLatticePointGenerator::operator()(Index ix) const {
  if (ix < 0 || ix >= this->m_total_lattice_points) {
    throw std::runtime_error(
        "Lattice point index out of range! Specified index " +
        std::to_string(ix) + " when there are " +
        std::to_string(m_total_lattice_points) + " lattice sites.");
  }

  auto mnp = this->_make_smith_normal_form_lattice_point(ix);
  return this->_normalize_lattice_point(mnp);
}

std::vector<UnitCell> make_lattice_points(
    const OrderedLatticePointGenerator::matrix_type &transformation_matrix) {
  std::vector<UnitCell> all_lattice_points;
  OrderedLatticePointGenerator generate_point(transformation_matrix);

  auto total_lattice_points = generate_point.size();

  for (int i = 0; i < total_lattice_points; ++i) {
    all_lattice_points.emplace_back(generate_point(i));
  }

  return all_lattice_points;
}

}  // namespace impl

//********************************************************************************************************************************//

std::vector<UnitCell> make_lattice_points(
    const Eigen::Matrix3l &transformation_matrix) {
  return impl::make_lattice_points(transformation_matrix);
}

std::vector<UnitCell> make_lattice_points(const Lattice &tiling_unit,
                                          const Lattice &superlattice,
                                          double tol) {
  Eigen::Matrix3l transformation_matrix =
      make_transformation_matrix_to_super(tiling_unit, superlattice, tol);
  return make_lattice_points(transformation_matrix);
}

}  // namespace xtal
}  // namespace CASM
