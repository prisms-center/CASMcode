#include "casm/crystallography/LatticePointWithin.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/external/Eigen/Core"
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

namespace CASM {
  namespace xtal {
    void LatticePointWithin::_throw_if_bad_transformation_matrix(const matrix_type &transformation_matrix) {
      if(transformation_matrix.determinant() == 0) {
        throw std::runtime_error(
          "The transformation matrix that converts the tiling unit to the superlattice is singular, and therefore not valid.");
      }
      return;
    }

    LatticePointWithin::matrix_type
    LatticePointWithin::_make_transformation_matrix(const Lattice &tiling_unit, const Lattice &superlattice, double tol) {
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

    UnitCellCoord LatticePointWithin::operator()(const UnitCellCoord &bijk) const {
      return UnitCellCoord(bijk.sublattice(), this->operator()(bijk.unitcell()));
    }

    //********************************************************************************************************************************//

    PrimGrid__::PrimGrid__(const matrix_type &transformation_matrix)
      : m_bring_within_f(transformation_matrix), m_total_lattice_points(transformation_matrix.determinant()) {
      smith_normal_form(transformation_matrix, this->m_smith_normal_U, this->m_smith_normal_S, this->m_smith_normal_V);

      auto S_diagonal = this->m_smith_normal_S.diagonal();

      m_stride[0] = S_diagonal[0];
      m_stride[1] = S_diagonal[0] * S_diagonal[1];

      assert(transformation_matrix == m_smith_normal_U * m_smith_normal_S * m_smith_normal_V);
    }

    PrimGrid__::vector_type PrimGrid__::_normalize_lattice_point(const vector_type &mnp) const {
      vector_type ijk(0, 0, 0);
      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          ijk[i] += m_smith_normal_U(i, j) * mnp[j];
        }
      }
      return m_bring_within_f(ijk);
    }

    PrimGrid__::vector_type PrimGrid__::_make_smith_normal_form_lattice_point(Index ix) const {
      vector_type mnp((ix % m_stride[1]) % m_stride[0], (ix % m_stride[1]) / m_stride[0], ix / m_stride[1]);
      return mnp;
    }

    PrimGrid__::vector_type PrimGrid__::operator()(Index ix) const {
      if(ix < 0 || ix >= this->m_total_lattice_points) {
        throw std::runtime_error("PrimGrid index out of range! Specified index " + std::to_string(ix) + " when there are " +
                                 std::to_string(m_total_lattice_points) + " lattice sites.");
      }

      auto mnp = this->_make_smith_normal_form_lattice_point(ix);
      return this->_normalize_lattice_point(mnp);
    }

  } // namespace xtal
} // namespace CASM
