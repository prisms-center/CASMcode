#ifndef SYMBASISPERMUTE_HH
#define SYMBASISPERMUTE_HH

#include "casm/container/Array.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/symmetry/SymOpRepresentation.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/casm_io/jsonParser.hh"

namespace CASM {
  class MasterSymGroup;
  class Lattice;

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
  ///     UnitCellSite(UnitCell(0,0,0), b) -> UnitCellSite'
  ///   which is used to set the sublat and is added to u' along with L.inv*T
  class SymBasisPermute: public SymOpRepresentation {
  private:

    /// Array of UnitCellCoords, of length basis.size(). Site (b,0,0,0) goes to
    /// (b_new,i,j,k)==m_ucc_permute[b] after application of symmetry.
    std::vector<UnitCellCoord> m_ucc_permute;

    /// Transform fractional coordinates, integer version of SymOp::matrix()
    Eigen::Matrix3l m_point_mat;


  public:

    typedef Index size_type;

    /// Construct SymBasisPermute
    template<typename StrucType>
    SymBasisPermute(const SymOp &op, const StrucType &struc, double tol);

    /// \brief Apply to a UnitCellCoord, in place
    UnitCellCoord &apply(UnitCellCoord &value) const {
      value.unitcell() = m_point_mat * value.unitcell() +
                         m_ucc_permute[value.sublat()].unitcell();
      value.sublat() = m_ucc_permute[value.sublat()].sublat();
      return value;
    }

    /// \brief Copy UnitCellCoord and apply
    UnitCellCoord copy_apply(UnitCellCoord value) const {
      return apply(value);
    }

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
    SymBasisPermute const *get_ucc_permutation() const override {
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

  /// \brief Apply symmetry to a UnitCellCoord
  template<typename BasisPermutable>
  UnitCellCoord &apply(const SymOp &op, UnitCellCoord &value, const BasisPermutable &obj);

  /// \brief Apply symmetry to a UnitCellCoord
  inline UnitCellCoord &apply(const SymBasisPermute &op, UnitCellCoord &value) {
    return op.apply(value);
  }


  // ---- SymBasisPermute Definitions --------------------

  /// Construct SymBasisPermute
  template<typename StrucType>
  SymBasisPermute::SymBasisPermute(const SymOp &op, const StrucType &struc, double tol) {
    SymOp::matrix_type frac_op(cart2frac(op.matrix(), struc.lattice()));
    if(!is_integer(frac_op, tol)) {
      throw std::runtime_error(
        std::string("Error in 'SymBasisPermute(const SymOp& op, const StrucType& struc, double tol)'\n") +
        "  Could not get integer point transformation matrix.");
    }

    m_point_mat = lround(frac_op);

    // Determine how basis sites transform from the origin unit cell
    for(int b = 0; b < struc.basis.size(); b++) {
      m_ucc_permute.push_back(UnitCellCoord(CASM::copy_apply(op, struc.basis[b]), struc, tol));
    }
  }


  /// \brief Apply symmetry to a UnitCellCoord
  ///
  /// \param op The symmetry operation to apply
  /// \param value The UnitCellCoord being transformed
  /// \param obj The object that the UnitCellCoord coordinates refer to, typically a primitive Structure
  ///
  /// - Requires BasisPermutable::basis_permutation_symrep_ID() to obtain the
  ///   SymBasisPermute representation
  template<typename BasisPermutable>
  UnitCellCoord &apply(const SymOp &op, UnitCellCoord &value, const BasisPermutable &obj) {

    return op.get_basis_permute_rep(obj.basis_permutation_symrep_ID())->apply(value);
  }

  /** @} */
}
#endif
