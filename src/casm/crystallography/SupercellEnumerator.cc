#include "casm/crystallography/SupercellEnumerator.hh"

#include <boost/math/special_functions/round.hpp>
#include "casm/external/Eigen/Dense"

#include "casm/crystallography/Structure.hh"

namespace CASM {

  template<>
  SupercellEnumerator<Lattice>::SupercellEnumerator(Lattice unit,
                                                    double tol,
                                                    size_type begin_volume,
                                                    size_type end_volume) :
    m_unit(unit),
    m_lat(unit),
    m_begin_volume(begin_volume),
    m_end_volume(end_volume) {

    m_lat.generate_point_group(m_point_group, tol);

  }

  template<>
  SupercellEnumerator<Lattice>::SupercellEnumerator(Lattice unit,
                                                    const SymGroup &point_grp,
                                                    size_type begin_volume,
                                                    size_type end_volume) :
    m_unit(unit),
    m_lat(unit),
    m_point_group(point_grp),
    m_begin_volume(begin_volume),
    m_end_volume(end_volume) {}

  /// \brief Return canonical hermite normal form of the supercell matrix, and op used to find it
  ///
  /// \returns std::pair<Eigen::Matrix3i, Eigen::Matrix3d> of H in canonical form, and op used to find it from T
  ///
  /// \param T a supercell matrix (Eigen::Matrix3i), such that S = U*T,
  ///          where S is the superlattice and U the unit lattice, as column vector matrices
  /// \param unitcell the unit BasicStructure<Site>
  ///
  /// Canonical form is such that T is in hermite normal form (as from casm::hermite_normal_form),
  ///  and the unrolled coefficients
  /// \code
  /// [a f e]
  /// [0 b d] -> abcdef
  /// [0 0 c]
  /// \endcode
  /// form the highest lexicographic order when considering equivalent superlattices by point group operations.
  ///
  /// - Equivalent superlattices can be obtained using point group operations: S' = op*S = U*H*V,
  ///   where V is integer and has determinant +/- 1
  /// - Substituting S = U*T, we have op*U*T = U*H*V.
  /// - Or H*V = U.inverse*op*U*T, which is the hermite normal form of U.inverse*op*U*T
  /// - So T is canonical if it is in hermite normal form and for all point group operations
  ///   it has a higher lexicographic order than the resulting H
  ///
  ///
  /// \relatesalso Lattice
  ///
  std::pair<Eigen::MatrixXi, Eigen::MatrixXd> canonical_hnf(const Eigen::MatrixXi &T, const BasicStructure<Site> &unitcell) {

    Eigen::Matrix3d lat = unitcell.lattice().lat_column_mat();
    Structure unitstruc(unitcell);
    SymGroup pg = unitstruc.point_group();

    Eigen::Matrix3i H, H_init, H_canon;

    int i_canon = 0;

    // get T in hermite normal form
    H = hermite_normal_form(T).first;
    H_canon = H;
    H_init = H;

    for(int i = 0; i < pg.size(); i++) {

      Eigen::Matrix3i transformed = iround(lat.inverse() * (Eigen::Matrix3d(pg[i].get_matrix(CART)) * lat)) * H_init;

      H = hermite_normal_form(transformed).first;

      // canonical only if H_canon is '>=' H, for all H, so if H '>' m_current, make H the H_canon
      if(H(0, 0) > H_canon(0, 0)) {
        i_canon = i;
        H_canon = H;
        continue;
      }
      if(H(0, 0) < H_canon(0, 0))
        continue;

      if(H(1, 1) > H_canon(1, 1)) {
        i_canon = i;
        H_canon = H;
        continue;
      }
      if(H(1, 1) < H_canon(1, 1))
        continue;

      if(H(2, 2) > H_canon(2, 2)) {
        i_canon = i;
        H_canon = H;
        continue;
      }
      if(H(2, 2) < H_canon(2, 2))
        continue;

      if(H(1, 2) > H_canon(1, 2)) {
        i_canon = i;
        H_canon = H;
        continue;
      }
      if(H(1, 2) < H_canon(1, 2))
        continue;

      if(H(0, 2) > H_canon(0, 2)) {
        i_canon = i;
        H_canon = H;
        continue;
      }
      if(H(0, 2) < H_canon(0, 2))
        continue;

      if(H(0, 1) > H_canon(0, 1)) {
        i_canon = i;
        H_canon = H;
        continue;
      }
      if(H(0, 1) < H_canon(0, 1))
        continue;

    }

    Eigen::Matrix3d op_canon = pg[i_canon].get_matrix(CART);
    return std::make_pair<Eigen::MatrixXi, Eigen::MatrixXd>(H_canon, op_canon);
  }

}
