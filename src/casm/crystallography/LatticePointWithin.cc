#include "casm/crystallography/LatticePointWithin.hh"
#include "casm/external/Eigen/Core"
#include <exception>
#include <stdexcept>

namespace CASM {
  namespace xtal {
    void LatticePointWithin::_throw_if_bad_transformation_matrix(const matrix_type &transformation_matrix) {
      if(transformation_matrix.determinant() == 0) {
        throw std::runtime_error(
          "The transformation matrix that converts the tiling unit to the superlattice is singular, and therefore not valid.");
      }

      return;
    }

    LatticePointWithin::matrix_type LatticePointWithin::_make_transformation_matrix(const Lattice &tiling_unit, const Lattice &superlattice, double tol) {
      Eigen::Matrix3d direct_transformation_matrix = tiling_unit.lat_column_mat().inverse() * superlattice.lat_column_mat();

      matrix_type rounded_transformation_matrix = round(direct_transformation_matrix).cast<long>();
      LatticePointWithin::_throw_if_bad_transformation_matrix(rounded_transformation_matrix);

      Eigen::Matrix3d matrix_error = direct_transformation_matrix - rounded_transformation_matrix.cast<double>();
      if(!matrix_error.isZero(tol)) {
        throw std::runtime_error("The provided tiling unit and superlattice are not related by an integer transformation.");
      }

      return rounded_transformation_matrix;
    }

    LatticePointWithin::vector_type LatticePointWithin::operator()(const Eigen::Vector3l &ijk) const {
      vector_type vec2 = this->m_transformation_matrix_adjugate * ijk;
      auto vol = this->m_total_lattice_points_in_superlattice;
      vec2[0] = ((vec2[0] % vol) + vol) % vol;
      vec2[1] = ((vec2[1] % vol) + vol) % vol;
      vec2[2] = ((vec2[2] % vol) + vol) % vol;
      return (this->m_transformation_matrix * vec2) / vol;
    }

  } // namespace xtal
} // namespace CASM
