#include "casm/crystallography/LatticePointWithin.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/external/Eigen/Core"
#include "casm/global/eigen.hh"
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

namespace CASM {
  namespace xtal {
    void LatticePointWithin_f::_throw_if_bad_transformation_matrix(const matrix_type &transformation_matrix) {
      if(transformation_matrix.determinant() == 0) {
        throw std::runtime_error(
          "The transformation matrix that converts the tiling unit to the superlattice is singular, and therefore not valid.");
      }
      return;
    }

    LatticePointWithin_f::vector_type LatticePointWithin_f::operator()(const vector_type &ijk) const {
      vector_type vec2 = this->m_transformation_matrix_adjugate * ijk;
      auto vol = this->m_total_lattice_points_in_superlattice;
      vec2[0] = ((vec2[0] % vol) + vol) % vol;
      vec2[1] = ((vec2[1] % vol) + vol) % vol;
      vec2[2] = ((vec2[2] % vol) + vol) % vol;
      return (this->m_transformation_matrix * vec2) / vol;
    }

    UnitCellCoord LatticePointWithin_f::operator()(const UnitCellCoord &bijk) const {
      return UnitCellCoord(bijk.sublattice(), this->operator()(bijk.unitcell()));
    }

    //********************************************************************************************************************************//

    OrderedLatticePointGenerator::OrderedLatticePointGenerator(const matrix_type &transformation_matrix)
      : m_bring_within_f(transformation_matrix), m_total_lattice_points(std::abs(transformation_matrix.determinant())) {
      smith_normal_form(transformation_matrix, this->m_smith_normal_U, this->m_smith_normal_S, this->m_smith_normal_V);

      auto S_diagonal = this->m_smith_normal_S.diagonal();

      m_stride[0] = S_diagonal[0];
      m_stride[1] = S_diagonal[0] * S_diagonal[1];

      assert(transformation_matrix == m_smith_normal_U * m_smith_normal_S * m_smith_normal_V);
    }

    OrderedLatticePointGenerator::vector_type OrderedLatticePointGenerator::_normalize_lattice_point(const vector_type &mnp) const {
      vector_type ijk(0, 0, 0);
      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          ijk[i] += m_smith_normal_U(i, j) * mnp[j];
        }
      }
      return m_bring_within_f(ijk);
    }

    OrderedLatticePointGenerator::vector_type OrderedLatticePointGenerator::_make_smith_normal_form_lattice_point(Index ix) const {
      vector_type mnp((ix % m_stride[1]) % m_stride[0], (ix % m_stride[1]) / m_stride[0], ix / m_stride[1]);
      return mnp;
    }

    OrderedLatticePointGenerator::vector_type OrderedLatticePointGenerator::operator()(Index ix) const {
      if(ix < 0 || ix >= this->m_total_lattice_points) {
        throw std::runtime_error("PrimGrid index out of range! Specified index " + std::to_string(ix) + " when there are " +
                                 std::to_string(m_total_lattice_points) + " lattice sites.");
      }

      auto mnp = this->_make_smith_normal_form_lattice_point(ix);
      return this->_normalize_lattice_point(mnp);
    }

    //********************************************************************************************************************************//

    std::vector<UnitCell> make_lattice_points(const OrderedLatticePointGenerator::matrix_type &transformation_matrix) {
      std::vector<UnitCell> all_lattice_points;
      OrderedLatticePointGenerator generate_point(transformation_matrix);

      auto total_lattice_points = generate_point.size();

      for(int i = 0; i < total_lattice_points; ++i) {
        all_lattice_points.emplace_back(std::move(generate_point(i)));
      }

      return all_lattice_points;
    }

    std::vector<UnitCell> make_lattice_points(const Eigen::Matrix3i &transformation_matrix) {
      return make_lattice_points(Eigen::Matrix3l(transformation_matrix.cast<long>()));
    }

    std::vector<UnitCell> make_lattice_points(const Lattice &tiling_unit, const Lattice &superlattice) {
      auto transformation_matrix = make_transformation_matrix(tiling_unit, superlattice, TOL);
      return make_lattice_points(transformation_matrix);
    }

  } // namespace xtal
} // namespace CASM
