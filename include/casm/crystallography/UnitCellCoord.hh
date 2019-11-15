#ifndef UNITCELLCOORD_HH
#define UNITCELLCOORD_HH

#include <iostream>

#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"
#include "casm/misc/Comparisons.hh"

namespace CASM {
  class jsonParser;
  class SymOp;

  namespace xtal {
    class Coordinate;
    class Site;
    template <typename CoordType>
    class BasicStructure;
    class Structure;
    class Lattice;
    class UnitCellCoord;
  }

  namespace sym {
    xtal::UnitCellCoord &apply(const CASM::SymOp &op, xtal::UnitCellCoord &ucc, const xtal::Structure &prim);
  }

  namespace xtal {

    /** \ingroup Coordinate
     *  @{
     */

    /// \brief Unit Cell Indices
    ///
    /// - Integer vector to represent a particular unit cell using multiples of the unit cell vectors
    /// - Is a <a href="http://eigen.tuxfamily.org/dox/group__QuickRefPage.html">Eigen::Vector3i</a>.
    ///
    class UnitCell : public Eigen::Vector3l {
    public:
      UnitCell(void) : Eigen::Vector3l() {}

      UnitCell(Index a, Index b, Index c) : Eigen::Vector3l(a, b, c) {}

      // This constructor allows you to construct MyVectorType from Eigen expressions
      template <typename OtherDerived>
      UnitCell(const Eigen::MatrixBase<OtherDerived> &other) : Eigen::Vector3l(other) {
      }

      // This method allows you to assign Eigen expressions to MyVectorType
      template <typename OtherDerived>
      UnitCell &operator=(const Eigen::MatrixBase<OtherDerived> &other) {
        this->Eigen::Vector3l::operator=(other);
        return *this;
      }

      /// \brief Compare UnitCell
      bool operator<(const UnitCell &B) const {
        const auto &A = *this;
        for(Index i = 0; i < 3; i++) {
          if(A(i) < B(i)) {
            return true;
          }
          if(A(i) > B(i)) {
            return false;
          }
        }
        return false;
      }
    };

    /// Convert lattice point a unitcell
    UnitCell make_unitcell(Coordinate const &lattice_point);

    /// \brief CRTP class to implement '-=', '+', and '-' in terms of '+='
    ///
    /// Requires:
    /// - MostDerived& MostDerived::operator+=(UnitCell frac)
    ///
    template <typename Base>
    struct Translatable : public Base {

      typedef typename Base::MostDerived MostDerived;
      using Base::derived;

      MostDerived &operator-=(UnitCell frac) {
        return derived() += -frac;
      }

      MostDerived operator+(UnitCell frac) const {
        MostDerived tmp{derived()};
        return tmp += frac;
      }

      MostDerived operator-(UnitCell frac) const {
        MostDerived tmp{derived()};
        return tmp += -frac;
      }
    };

    /* -- UnitCellCoord Declarations ------------------------------------- */

    /// \brief Unit Cell Coordinates
    ///
    /// - Represent a crystal site using UnitCell indices and sublattice index
    ///
    class UnitCellCoord : public Comparisons<Translatable<CRTPBase<UnitCellCoord>>> {

    public:
      typedef BasicStructure<Site> PrimType;

      UnitCellCoord(Index _sublat, const UnitCell &_unitcell) : m_unitcell(_unitcell), m_sublat(_sublat) {
        if(this->sublattice() < 0) {
          throw std::runtime_error("Error in UnitCellCoord. Construction requires a positive sublattice index.");
        }
      }

      UnitCellCoord(Index _sublat, Index i, Index j, Index k) : UnitCellCoord(_sublat, UnitCell(i, j, k)) {}

      explicit UnitCellCoord() : UnitCellCoord(0, 0, 0, 0) {}

      UnitCellCoord(const Coordinate &coord, double tol);

      UnitCellCoord(const UnitCellCoord &B) = default;

      UnitCellCoord &operator=(const UnitCellCoord &B) = default;

      UnitCellCoord(UnitCellCoord &&B) = default;

      UnitCellCoord &operator=(UnitCellCoord &&B) = default;

      const UnitCell &unitcell() const {
        return m_unitcell;
      }

      Index unitcell(Index i) const {
        return m_unitcell[i];
      }

      Index sublattice() const {
        return m_sublat;
      }

      //TODO: This should take a Lattice not a PrimType
      /// \brief Get corresponding coordinate
      Coordinate coordinate(const PrimType &prim) const;

      /// \brief Get a copy of corresponding site
      Site site(const PrimType &prim) const;

      /// \brief Get reference to corresponding sublattice site in the unit structure
      const Site &sublattice_site(const PrimType &prim) const;

      Index operator[](Index i) const;

      UnitCellCoord &operator+=(UnitCell frac);

      bool operator<(const UnitCellCoord &B) const;

      /* UnitCellCoord& apply_sym(const CASM::SymOp& op); */

      /* UnitCellCoord copy_apply(const CASM::SymOp& op) const; */

    private:
      /// make _eq accessible
      friend struct Comparisons<Translatable<CRTPBase<UnitCellCoord>>>;
      //Grant access to in place applying symmetry
      friend UnitCellCoord &sym::apply(const CASM::SymOp &op, UnitCellCoord &ucc, const Structure &prim);

      UnitCell &_unitcell() {
        return m_unitcell;
      }

      Index &_sublattice() {
        return m_sublat;
      }

      bool eq_impl(const UnitCellCoord &B) const;

      /// Returns false if the sublattice of *this is greater than the number of sites in the given prmitive structure
      bool _is_compatible_with_prim(const PrimType &prim) const;

      // TODO: Should this be made into an actual exception class?
      static void _throw_incompatible_primitive_cell();

      UnitCell m_unitcell;
      Index m_sublat;
    };

    inline std::ostream &operator<<(std::ostream &sout, const UnitCellCoord &site) {
      return sout << site.sublattice() << ", " << site.unitcell().transpose();
    }

    inline Index UnitCellCoord::operator[](Index i) const {
      if(i == 0) {
        return m_sublat;
      }
      return unitcell(i - 1);
    }

    inline UnitCellCoord &UnitCellCoord::operator+=(UnitCell frac) {
      m_unitcell += frac;
      return *this;
    }

    /// \brief Compare UnitCellCoord
    inline bool UnitCellCoord::operator<(const UnitCellCoord &B) const {
      const auto &A = *this;
      for(Index i = 0; i < 3; i++) {
        if(A.unitcell()(i) < B.unitcell()(i)) {
          return true;
        }
        if(A.unitcell()(i) > B.unitcell()(i)) {
          return false;
        }
      }
      if(A.sublattice() < B.sublattice()) {
        return true;
      }

      return false;
    }

    inline bool UnitCellCoord::eq_impl(const UnitCellCoord &B) const {
      const auto &A = *this;
      return A.unitcell()(0) == B.unitcell()(0) && A.unitcell()(1) == B.unitcell()(1) &&
             A.unitcell()(2) == B.unitcell()(2) && A.sublattice() == B.sublattice();
    }

    /** @} */
  } // namespace xtal

  /// \brief Print to json as [b, i, j, k]
  jsonParser &to_json(const xtal::UnitCellCoord &ucc_val, jsonParser &fill_json);

  /// \brief Read from json [b, i, j, k]
  void from_json(xtal::UnitCellCoord &fill_value, const jsonParser &read_json);

} // namespace CASM
#endif

