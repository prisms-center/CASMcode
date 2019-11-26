#include "casm/crystallography/UnitCellWithin.hh"
#include "casm/external/Eigen/src/Core/Matrix.h"
#include <exception>
#include <stdexcept>

namespace CASM {
  namespace xtal {
    void UnitCellWithin::_throw_if_bad_transformation_matrix(const matrix_type &transformation_matrix) {
      if(transformation_matrix.determinant() == 0) {
        throw std::runtime_error(
          "The transformation matrix that converts the tiling unit to the superlattice is singular, and therefore not valid.");
      }

      return;
    }

    UnitCellWithin::matrix_type UnitCellWithin::_make_transformation_matrix(const Lattice &tiling_unit, const Lattice &superlattice, double tol) {
      Eigen::Matrix3d direct_transformation_matrix = tiling_unit.lat_column_mat().inverse() * superlattice.lat_column_mat();

      matrix_type rounded_transformation_matrix = round(direct_transformation_matrix).cast<long>();
      UnitCellWithin::_throw_if_bad_transformation_matrix(rounded_transformation_matrix);

      Eigen::Matrix3d matrix_error = direct_transformation_matrix - rounded_transformation_matrix.cast<double>();
      if(!matrix_error.isZero(tol)) {
        throw std::runtime_error("The provided tiling unit and superlattice are not related by an integer transformation.");
      }

      return rounded_transformation_matrix;
    }


    /* { */
    /*       vector_type vec2 = m_plane_mat * ijk; */

    /*       vec2[0] = ((vec2[0] % m_N_vol) + m_N_vol) % m_N_vol; */
    /*       vec2[1] = ((vec2[1] % m_N_vol) + m_N_vol) % m_N_vol; */
    /*       vec2[2] = ((vec2[2] % m_N_vol) + m_N_vol) % m_N_vol; */

    /*       return (m_trans_mat * vec2) / m_N_vol; */
    /* } */

  } // namespace xtal
} // namespace CASM
