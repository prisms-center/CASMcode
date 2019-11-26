#ifndef UNITCELLWITHIN_HH
#define UNITCELLWITHIN_HH

#include "casm/crystallography/Lattice.hh"
#include "casm/external/Eigen/Core"
#include "casm/external/Eigen/src/Core/Matrix.h"
#include "casm/global/eigen.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
  namespace xtal {
    class UnitCell;

    /**
     * Handles bringing a UnitCell (i,j,k values) within a particular
     * superlattice. Provide the tiling unit and superlattice at construction,
     * and use this class to bring any external UnitCell within the
     * superlattice
     */

    struct UnitCellWithin {
      typedef Eigen::Matrix<long, 3, 3> matrix_type;
      typedef Eigen::Matrix<long, 3, 1> vector_type;

      /// Specify the tiling unit and a transformation that turns the tiling unit into the desired
      /// superlattice into which UnitCells should be brought within
      UnitCellWithin(const matrix_type &superlattice_transformation_matrix)
        : m_transformation_matrix(superlattice_transformation_matrix),
          m_transformation_matrix_adjugate(adjugate(superlattice_transformation_matrix)),
          m_total_lattice_points_in_superlattice(superlattice_transformation_matrix.determinant()) {
        _throw_if_bad_transformation_matrix(this->m_transformation_matrix);
      }

      UnitCellWithin(const Eigen::Matrix3i &superlattice_transformation_matrix):
        UnitCellWithin(matrix_type(superlattice_transformation_matrix.cast<long>()))
      {}

      /// Specify the tiling unit, and the superlattice into which UnitCells should be brought within.
      /// The superlattice must be an integer transformation of the tiling unit
      UnitCellWithin(const Lattice &tiling_unit, const Lattice &superlattice):
        UnitCellWithin(UnitCellWithin::_make_transformation_matrix(tiling_unit, superlattice, TOL))
      {}

      ///Brings the given lattice point within the superlattice
      Eigen::Vector3l operator()(const Eigen::Vector3l &ijk) const;
      UnitCell operator()(const UnitCell &ijk) const;

    private:
      /// Integer matrix that converts the tiling unit into the superlattice.
      /// For a tiling unit U, superlattice S, and transformation matrix T -> S=U*T
      matrix_type m_transformation_matrix;

      /// The adjugate matrix of the transformation matrix.
      /// adjugate(T)=det(T)*inverse(T)
      matrix_type m_transformation_matrix_adjugate;

      long int m_total_lattice_points_in_superlattice;

      /// Calculates the transformation matrix that takes the tiling unit to the superlattice.
      /// Throws exceptions if the superlattice isn't compatible with its tiling unit
      static matrix_type _make_transformation_matrix(const Lattice &tiling_unit, const Lattice &superlattice, double tol);

      /// Throws exception if the transformation matrix has determinant 0
      static void _throw_if_bad_transformation_matrix(const matrix_type &transformation_matrix);
    };
  } // namespace xtal
} // namespace CASM

#endif
