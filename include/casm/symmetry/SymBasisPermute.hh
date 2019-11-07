#ifndef SYMBASISPERMUTE_HH
#define SYMBASISPERMUTE_HH

#include <vector>
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/symmetry/SymOpRepresentation.hh"

namespace CASM {
  namespace xtal {
    class Lattice;
    class UnitCellCoord;
  }
  using xtal::Lattice;
  using xtal::UnitCellCoord;

  class MasterSymGroup;
  class jsonParser;

  /** \ingroup SymOp
   *  @{
   */

  ///\brief  SymBasisPermute describes how a symmetry operation permutes atoms in a basis
  ///
  /// - ::SymOp transforms Cartesian coordinate (x) like:
  ///     x' = R*x + T, where R is a point transformation matrix, and T a translation vector
  /// - SymBasisPermute transforms UnitCell (u) without basis like:
  ///    Lu' = R*L*u + T
  ///     u' = L.inv*R*L*u + L.inv*T, where L is the lattice, as a column vector matrix
  /// - For transforming basis sites, a lookup table is stored that maps
  ///     UnitCellCoord(struc, b, UnitCell(0,0,0)) -> UnitCellCoord'
  ///   which is used to set the sublat and is added to u' along with L.inv*T
  /// - L.inv*R*L is stored as SymBasisPermute::matrix()
  /// - The result of transforming (b,0,0,0) is stored in SymBasisPermute::data()[b]
  /// - So a UnitCellCoord, u, is transformed by a SymOp, op, like:
  ///   \code
  ///   const SymBasisPermute &rep = *op.get_basis_permute_rep(u.unit().basis_permutation_symrep_ID());
  ///   u.unitcell() = rep.matrix()*u.unitcell() + rep.data()[u.sublat()].unitcell();
  ///   u.sublat() = rep.data()[u.sublat()].sublat();
  ///   \endcode
  /// - Additional translations, such as those in supercell factor group operations,
  ///   may be stored in SymOp::integral_tau(). Those additional contributions
  ///   can be included with:
  ///   - u.unitcell() += (u.unit().lattice().inv_lat_column_mat() * op.integral_tau()).cast<long>();
  ///
  class SymBasisPermute: public SymOpRepresentation {
  private:

    /// vector of UnitCellCoords, of length basis.size(). Site (b,0,0,0) goes to
    /// (b_new,i,j,k)==m_ucc_permute[b] after application of symmetry.
    std::vector<UnitCellCoord> m_ucc_permute;

    /// Transform fractional coordinates, integer version of SymOp::matrix()
    Eigen::Matrix3l m_point_mat;


  public:

    typedef Index size_type;

    /// Construct SymBasisPermute
    SymBasisPermute(SymOp const &_op, Lattice const &_lat,  std::vector<UnitCellCoord> const &_ucc_permute);

    /// Return pointer to a copy of this SymBasisPermute
    SymOpRepresentation *copy() const override {
      return new SymBasisPermute(*this);
    };

    /// Return number of basis sites
    size_type size() const {
      return m_ucc_permute.size();
    }

    /// Return UnitCellCoord that (b,0,0,0) transforms to
    UnitCellCoord operator[](Index b) const {
      return m_ucc_permute[b];
    }

    /// Return UnitCellCoord that (b,0,0,0) transforms to
    UnitCellCoord at(Index b) const {
      return m_ucc_permute.at(b);
    }

    /// Get this from a SymOp
    SymBasisPermute const *ucc_permutation() const override {
      return this;
    };

    /// Get underlying data (data()[b] is the result of transforming (b,0,0,0))
    const std::vector<UnitCellCoord> &data() const {
      return m_ucc_permute;
    }

    /// Get underlying integer transformation amtrix
    const Eigen::Matrix3l &matrix() const {
      return m_point_mat;
    }

    jsonParser &to_json(jsonParser &json) const override {
      return json;
    };
    void from_json(const jsonParser &json) override {};

  };

}
#endif
