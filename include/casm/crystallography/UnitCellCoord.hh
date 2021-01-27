#ifndef UNITCELLCOORD_HH
#define UNITCELLCOORD_HH

#include <iostream>
#include <stdexcept>

#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"
#include "casm/misc/Comparisons.hh"

namespace CASM {
namespace xtal {
class Coordinate;
class Site;
class BasicStructure;
class Lattice;
class UnitCellCoord;
class Superlattice;
}  // namespace xtal

namespace xtal {

/** \ingroup Coordinate
 *  @{
 */

/// \brief Unit Cell Indices
///
/// - Integer vector to represent a particular unit cell using multiples of the
/// unit cell vectors
/// - Is a <a
/// href="http://eigen.tuxfamily.org/dox/group__QuickRefPage.html">Eigen::Vector3i</a>.
///
class UnitCell : public Eigen::Vector3l {
 public:
  UnitCell(void) : Eigen::Vector3l() {}

  /// Construct from integral fractional values, relative to the tiling unit.
  template <typename OtherDerived>
  UnitCell(const Eigen::MatrixBase<OtherDerived> &integral_coordinate)
      : Eigen::Vector3l(integral_coordinate) {}

  UnitCell(Index a, Index b, Index c) : UnitCell(Eigen::Vector3l(a, b, c)) {}

  /// Convert lattice point to a unitcell
  static UnitCell from_coordinate(Coordinate const &lattice_point);

  /// Convert a Cartesian coordinate into a unitcell by rounding fractional
  /// coordinates to the provided lattice
  static UnitCell from_cartesian(const Eigen::Vector3d &cartesian_coord,
                                 const Lattice &tiling_unit);

  /// Convert a unitcell to a lattice point coordinate, given the primitive
  /// tiling unit lattice
  Coordinate coordinate(const Lattice &tiling_unit) const;

  /// Finds a new UnitCell with values relative to the given tiling unit
  UnitCell reset_tiling_unit(const Lattice &current_tiling_unit,
                             const Lattice &new_tiling_unit) const;

  // This method allows you to assign Eigen expressions to MyVectorType
  template <typename OtherDerived>
  UnitCell &operator=(const Eigen::MatrixBase<OtherDerived> &other) {
    this->Eigen::Vector3l::operator=(other);
    return *this;
  }

  /// \brief Compare UnitCell
  bool operator<(const UnitCell &B) const {
    const auto &A = *this;
    for (Index i = 0; i < 3; i++) {
      if (A(i) < B(i)) {
        return true;
      }
      if (A(i) > B(i)) {
        return false;
      }
    }
    return false;
  }

 private:
  /// Throws exception to indicate that finding integral values resulted in
  /// errors much larger than the relevant lattice tolerance
  static void _throw_large_rounding_error() {
    throw std::runtime_error(
        "Could not round values to integers within a reasonable tolerance");
  }
};

/// \brief CRTP class to implement '-=', '+', and '-' in terms of '+='
///
/// Requires:
/// - MostDerived& MostDerived::operator+=(UnitCell frac)
///
template <typename Base>
struct Translatable : public Base {
  typedef typename Base::MostDerived MostDerived;
  using Base::derived;

  MostDerived &operator-=(UnitCell frac) { return derived() += -frac; }

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
class UnitCellCoord
    : public Comparisons<Translatable<CRTPBase<UnitCellCoord>>> {
 public:
  typedef BasicStructure PrimType;

  UnitCellCoord(Index _sublat, const UnitCell &_unitcell)
      : m_unitcell(_unitcell), m_sublat(_sublat) {
    if (!valid_index(_sublat)) {
      throw std::runtime_error(
          "Error in UnitCellCoord. Construction requires a positive sublattice "
          "index.");
    }
  }

  UnitCellCoord(Index _sublat, Index i, Index j, Index k)
      : UnitCellCoord(_sublat, UnitCell(i, j, k)) {}

  explicit UnitCellCoord() : UnitCellCoord(0, 0, 0, 0) {}

  static UnitCellCoord from_coordinate(const PrimType &,
                                       const Coordinate &coord, double tol);

  UnitCellCoord(const UnitCellCoord &B) = default;

  UnitCellCoord &operator=(const UnitCellCoord &B) = default;

  UnitCellCoord(UnitCellCoord &&B) = default;

  UnitCellCoord &operator=(UnitCellCoord &&B) = default;

  const UnitCell &unitcell() const { return m_unitcell; }

  Index sublattice() const { return m_sublat; }

  /// \brief Get corresponding coordinate
  Coordinate coordinate(const PrimType &prim) const;

  /// \brief Get a copy of corresponding site
  Site site(const PrimType &prim) const;

  /// \brief Get reference to corresponding sublattice site in the unit
  /// structure
  const Site &sublattice_site(const PrimType &prim) const;

  Index operator[](Index i) const;

  UnitCellCoord &operator+=(UnitCell frac);

  bool operator<(const UnitCellCoord &B) const;

 private:
  /// make _eq accessible
  friend struct Comparisons<Translatable<CRTPBase<UnitCellCoord>>>;

  UnitCell &_unitcell() { return m_unitcell; }

  Index &_sublattice() { return m_sublat; }

  bool eq_impl(const UnitCellCoord &B) const;

  /// Returns false if the sublattice of *this is greater than the number of
  /// sites in the given prmitive structure
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
  if (i == 0) {
    return m_sublat;
  }
  return unitcell()[i - 1];
}

inline UnitCellCoord &UnitCellCoord::operator+=(UnitCell frac) {
  m_unitcell += frac;
  return *this;
}

/// \brief Compare UnitCellCoord
inline bool UnitCellCoord::operator<(const UnitCellCoord &B) const {
  const auto &A = *this;
  for (Index i = 0; i < 3; i++) {
    if (A.unitcell()(i) < B.unitcell()(i)) {
      return true;
    }
    if (A.unitcell()(i) > B.unitcell()(i)) {
      return false;
    }
  }
  if (A.sublattice() < B.sublattice()) {
    return true;
  }

  return false;
}

inline bool UnitCellCoord::eq_impl(const UnitCellCoord &B) const {
  const auto &A = *this;
  return A.unitcell()(0) == B.unitcell()(0) &&
         A.unitcell()(1) == B.unitcell()(1) &&
         A.unitcell()(2) == B.unitcell()(2) && A.sublattice() == B.sublattice();
}

/** @} */

/// Converts the position of a lattice site into a Coordinate within the
/// superlattice
Coordinate make_superlattice_coordinate(const UnitCell &ijk,
                                        const Superlattice &superlattice);
Coordinate make_superlattice_coordinate(const UnitCell &ijk,
                                        const Lattice &tiling_unit,
                                        const Lattice &superlattice);

}  // namespace xtal
}  // namespace CASM

/* #include "casm/crystallography/Coordinate.hh" */
namespace CASM {
namespace xtal {}
}  // namespace CASM

namespace std {
template <>
struct hash<CASM::xtal::UnitCell> {
  std::size_t operator()(const CASM::xtal::UnitCell &value) const;
};

template <>
struct hash<CASM::xtal::UnitCellCoord> {
  std::size_t operator()(const CASM::xtal::UnitCellCoord &value) const;
};
}  // namespace std
#endif
