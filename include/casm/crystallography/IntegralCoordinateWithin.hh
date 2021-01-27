#ifndef UNITCELLWITHIN_HH
#define UNITCELLWITHIN_HH

#include <array>
#include <vector>

#include "casm/external/Eigen/Core"
#include "casm/global/eigen.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
namespace xtal {
class UnitCell;
class UnitCellCoord;

/**
 * Handles bringing an integral coordinate (i,j,k value) within a particular
 * superlattice. Provide the tiling unit and superlattice at construction,
 * and use this class to bring any external UnitCell or UnitCellCoord within the
 * superlattice.
 */

class IntegralCoordinateWithin_f {
 public:
  typedef Eigen::Matrix<long, 3, 3> matrix_type;
  typedef Eigen::Matrix<long, 3, 1> vector_type;

  /// Specify the integer transformation matrix that turns the tiling unit
  /// into the desired superlattice. Lattice points (UnitCell) with fractional
  /// coordinates relative to the tiling unit will then be brought into
  /// the superlattice.
  IntegralCoordinateWithin_f(
      const matrix_type &superlattice_transformation_matrix)
      : m_transformation_matrix(superlattice_transformation_matrix),
        m_transformation_matrix_adjugate(
            adjugate(superlattice_transformation_matrix)) {
    _throw_if_bad_transformation_matrix(this->m_transformation_matrix);
  }

  IntegralCoordinateWithin_f(
      const Eigen::Matrix3i &superlattice_transformation_matrix)
      : IntegralCoordinateWithin_f(
            matrix_type(superlattice_transformation_matrix.cast<long>())) {}

  /// Brings the given lattice point within the superlattice
  vector_type operator()(const vector_type &ijk) const;

  template <typename UnitCellType>
  UnitCellType operator()(const UnitCellType &ijk) const {
    return UnitCellType(this->operator()(static_cast<vector_type>(ijk)));
  }

  UnitCellCoord operator()(const UnitCellCoord &bijk) const;

 private:
  /// Integer matrix that converts the tiling unit into the superlattice.
  /// For a tiling unit U, superlattice S, and transformation matrix T -> S=U*T
  matrix_type m_transformation_matrix;

  /// The adjugate matrix of the transformation matrix.
  /// adjugate(T)=det(T)*inverse(T)
  matrix_type m_transformation_matrix_adjugate;

  /// Throws exception if the transformation matrix has determinant 0
  static void _throw_if_bad_transformation_matrix(
      const matrix_type &transformation_matrix);
};

//********************************************************************************************************************************//

namespace impl {

/**
 * The Smith Normal Form of a transformation matrix allows the enumeration of
 * lattice points to be generated in a particular order. This class can generate
 * any lattice point given its index in the list in constant time.
 *
 * If you want to convert between linear index and bijk (UnitCellCoord) quickly,
 * you probably want to be using the UnitCellCoordIndexConverter class instead
 * of this.
 *
 * The original algorithm for constant time evaluation of the linear index that
 * makes this class possible was developed by John C. Thomas and is described
 * below:
 *
 * The transformation matrix T given at construction is decomposed using Smith
 * Normal Form decomposition, such that
 *
 * T=U*S*V
 *
 * where S is a diagonal matrix, and T converts the primitive lattice P to the
 * superlattice L such that
 *
 * L=P*T
 *
 * Combining the two equations we find that
 *
 * L*V.inv()=P*U*S
 *
 * In other words, P*U is a primitive lattice that perfectly tiles the
 * equivalent superlattice L*V.inv(), with the transformation matrix being S. If
 * we define L'=L*V.inv(), P'=P*U and T'=S, the realationship looks just like a
 * normal superlattice transformation:
 *
 * L'=P'*T'
 *
 * Since S (=T') is diagonal, points (m,n,p) on the transformed grid are a
 * convenient choice for indexing. Defining S0=S(0,0), S1=S(1,1), and S2=S(2,2),
 * an index i results in the following (m,n,p) values:
 *
 * m=(i%(S0*S1))%S0
 * n=(i%(S0*S1))/S0
 * p=i/(S0*S1)
 *
 * Once the (m,n,p) value is determined, it can be transformed back to the
 * original grid, to yield the desired (i,j,k) value, which is in terms of the
 * original superlattice tiling unit:
 *
 * (i,j,k)=U*(m,n,p)
 */

class OrderedLatticePointGenerator {
 public:
  typedef IntegralCoordinateWithin_f::matrix_type matrix_type;
  typedef IntegralCoordinateWithin_f::vector_type vector_type;

  /// Construct with the transformation matrix that takes the tiling unit to the
  /// superlattice
  OrderedLatticePointGenerator(const matrix_type &transformation_matrix);

  /// Given the index into the list of lattice points, return its corresponding
  /// lattice point, with guaranteed order
  vector_type operator()(Index ix) const;

  /// Returns the total number of unique lattice points that can be generated
  long int size() const { return m_total_lattice_points; }

 private:
  long int m_total_lattice_points;

  // Can map ijk values within the supercell
  IntegralCoordinateWithin_f m_bring_within_f;

  /// The Smith Normal Form decomposition is: trans_mat = U*S*V, with
  /// det(U)=det(V)=1; S is diagonal
  matrix_type m_smith_normal_U;
  /// The Smith Normal Form decomposition is: trans_mat = U*S*V, with
  /// det(U)=det(V)=1; S is diagonal
  matrix_type m_smith_normal_S;
  /// The Smith Normal Form decomposition is: trans_mat = U*S*V, with
  /// det(U)=det(V)=1; S is diagonal
  matrix_type m_smith_normal_V;

  /// stride maps canonical 3d index (m,n,p) onto linear index -- l = m +
  /// n*stride[0] + p*stride[1]
  std::array<int, 2> m_stride;

  /// Convert UnitCell created with strides (mnp) to a "normal" UnitCell in
  /// terms of the superlattice (ijk) U*mnp = ijk
  vector_type _normalize_lattice_point(const vector_type &mnp) const;

  /// Create a lattice point from the provided index, in the Smith Normal Form
  /// diagonalized space
  vector_type _make_smith_normal_form_lattice_point(Index ix) const;
};

/// Returns all the lattice points that exists within the superlattice that has
/// the given transformation matrix
std::vector<UnitCell> make_lattice_points(
    const OrderedLatticePointGenerator::matrix_type &transformation_matrix);
}  // namespace impl

//********************************************************************************************************************************//

std::vector<UnitCell> make_lattice_points(
    const Eigen::Matrix3l &transformation_matrix);

class Lattice;
/// Returns all the lattice points that exists when tiling the tiling unit
/// inside the superlattice
std::vector<UnitCell> make_lattice_points(const Lattice &tiling_unit,
                                          const Lattice &superlattice,
                                          double tol);

}  // namespace xtal
}  // namespace CASM

#endif
