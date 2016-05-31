#ifndef UNITCELLCOORD_HH
#define UNITCELLCOORD_HH

#include <iostream>

#include "casm/external/Eigen/Dense"
#include "casm/CASM_global_definitions.hh"
#include "casm/container/LinearAlgebra.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"

namespace CASM {

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
    template<typename OtherDerived>
    UnitCell(const Eigen::MatrixBase<OtherDerived> &other) :
      Eigen::Vector3l(other) {}

    // This method allows you to assign Eigen expressions to MyVectorType
    template<typename OtherDerived>
    UnitCell &operator=(const Eigen::MatrixBase <OtherDerived> &other) {
      this->Eigen::Vector3l::operator=(other);
      return *this;
    }
  };

  /* -- UnitCellCoord Declarations ------------------------------------- */

  /// \brief Unit Cell Coordinates
  ///
  /// - Represent a crystal site using UnitCell indices and sublattice index
  ///
  class UnitCellCoord {

  public:

    typedef BasicStructure<Site> UnitType;

    UnitCellCoord(const UnitType &unit);

    UnitCellCoord(const UnitType &unit, Index _sublat, const UnitCell &_unitcell);

    UnitCellCoord(const UnitType &unit, Index _sublat, Index i, Index j, Index k);

    UnitCellCoord(const UnitType &unit, const Coordinate &coord, double tol);


    UnitCellCoord(const UnitCellCoord &B) = default;

    UnitCellCoord &operator=(const UnitCellCoord &B) = default;

    UnitCellCoord(UnitCellCoord &&B) = default;

    UnitCellCoord &operator=(UnitCellCoord && B) = default;


    /// \brief Get unit structure reference
    const UnitType &unit() const;

    /// \brief Change unit structure, keeping indices constant
    void set_unit(const UnitType &_unit);

    /// \brief Get corresponding coordinate
    Coordinate coordinate() const;

    /// \brief Get corresponding site
    Site site() const;

    UnitCell &unitcell();
    const UnitCell &unitcell() const;

    Index &unitcell(Index i);
    const Index &unitcell(Index i) const;

    Index &sublat();
    const Index &sublat() const;

    Index &operator[](Index i);
    const Index &operator[](Index i) const;

    UnitCellCoord &operator+=(UnitCell frac);

    UnitCellCoord &operator-=(UnitCell frac);


  private:

    const UnitType *m_unit;
    UnitCell m_unitcell;
    Index m_sublat;

  };

  inline std::ostream &operator<<(std::ostream &sout, const UnitCellCoord &site) {
    return sout << site.sublat() << ", " << site.unitcell().transpose();
  }

  /// \brief Add UnitCell to UnitCellCoord
  inline UnitCellCoord operator+(UnitCell frac, UnitCellCoord site);

  /// \brief Add UnitCell to UnitCellCoord
  inline UnitCellCoord operator+(UnitCellCoord site, UnitCell frac);

  /// \brief Subtract UnitCell from UnitCellCoord
  inline UnitCellCoord operator-(UnitCellCoord site, UnitCell frac);

  /// \brief Compare UnitCellCoord
  inline bool operator<(const UnitCellCoord &A, const UnitCellCoord &B);

  /// \brief Compare UnitCellCoord
  inline bool operator>(const UnitCellCoord &A, const UnitCellCoord &B);

  /// \brief Compare UnitCellCoord
  inline bool operator<=(const UnitCellCoord &A, const UnitCellCoord &B);

  /// \brief Compare UnitCellCoord
  inline bool operator>=(const UnitCellCoord &A, const UnitCellCoord &B);

  /// \brief Compare UnitCellCoord
  inline bool operator==(const UnitCellCoord &A, const UnitCellCoord &B);

  /// \brief Compare UnitCellCoord
  inline bool operator!=(const UnitCellCoord &A, const UnitCellCoord &B);

  /// \brief Print to json as [b, i, j, k]
  jsonParser &to_json(const UnitCellCoord &ucc_val, jsonParser &fill_json);

  /// \brief Read from json [b, i, j, k]
  void from_json(UnitCellCoord &fill_value, const jsonParser &read_json);


  /* -- UnitCellCoord Definitions ------------------------------------- */

  inline UnitCellCoord::UnitCellCoord(const UnitType &unit) :
    m_unit(&unit) {}

  inline UnitCellCoord::UnitCellCoord(const UnitType &unit, Index _sublat, const UnitCell &_unitcell) :
    m_unit(&unit),
    m_unitcell(_unitcell),
    m_sublat(_sublat) {}

  inline UnitCellCoord::UnitCellCoord(const UnitType &unit, Index _sublat, Index i, Index j, Index k) :
    m_unit(&unit),
    m_unitcell(i, j, k),
    m_sublat(_sublat) {}

  inline UnitCellCoord::UnitCellCoord(const UnitType &unit, const Coordinate &coord, double tol) :
    m_unit(&unit) {
    for(Index b = 0; b < unit.basis.size(); ++b) {
      auto diff = coord - unit.basis[b];
      if(is_integer(diff.const_frac(), tol)) {
        *this = UnitCellCoord(unit, b, lround(diff.const_frac()));
        return;
      }
    }

    throw std::runtime_error(
      "Error in 'UnitCellCoord(CoordType coord, const StrucType& struc, double tol)'\n"
      "  No matching basis site found.");
  }

  inline const UnitCellCoord::UnitType &UnitCellCoord::unit() const {
    return *m_unit;
  }

  /// \brief Change unit structure, keeping indices constant
  inline void UnitCellCoord::set_unit(const UnitType &_unit) {
    m_unit = &_unit;
  }

  /// \brief Get corresponding coordinate
  inline Coordinate UnitCellCoord::coordinate() const {
    return unit().get_site(*this);
  }

  /// \brief Get corresponding site
  inline Site UnitCellCoord::site() const {
    return unit().get_site(*this);
  }

  inline UnitCell &UnitCellCoord::unitcell() {
    return m_unitcell;
  }

  inline const UnitCell &UnitCellCoord::unitcell() const {
    return m_unitcell;
  }

  inline Index &UnitCellCoord::unitcell(Index i) {
    return m_unitcell[i];
  }

  inline const Index &UnitCellCoord::unitcell(Index i) const {
    return m_unitcell[i];
  }

  inline Index &UnitCellCoord::sublat() {
    return m_sublat;
  }

  inline const Index &UnitCellCoord::sublat() const {
    return m_sublat;
  }

  inline Index &UnitCellCoord::operator[](Index i) {
    if(i == 0) {
      return m_sublat;
    }
    return unitcell(i - 1);
  }

  inline const Index &UnitCellCoord::operator[](Index i) const {
    if(i == 0) {
      return m_sublat;
    }
    return unitcell(i - 1);
  }

  inline UnitCellCoord &UnitCellCoord::operator+=(UnitCell frac) {
    m_unitcell += frac;
    return *this;
  }

  inline UnitCellCoord &UnitCellCoord::operator-=(UnitCell frac) {
    m_unitcell -= frac;
    return *this;
  }

  /// \brief Add UnitCell to UnitCellCoord
  inline UnitCellCoord operator+(UnitCell frac, UnitCellCoord site) {
    return site += frac;
  }

  /// \brief Add UnitCell to UnitCellCoord
  inline UnitCellCoord operator+(UnitCellCoord site, UnitCell frac) {
    return site += frac;
  }

  /// \brief Subtract UnitCell from UnitCellCoord
  inline UnitCellCoord operator-(UnitCellCoord site, UnitCell frac) {
    return site -= frac;
  }

  /// \brief Compare UnitCellCoord
  inline bool operator<(const UnitCellCoord &A, const UnitCellCoord &B) {
    for(Index i = 0; i < 3; i++) {
      if(A.unitcell()(i) < B.unitcell()(i)) {
        return true;
      }
      if(A.unitcell()(i) > B.unitcell()(i)) {
        return false;
      }
    }
    if(A.sublat() < B.sublat()) {
      return true;
    }

    return false;
  }

  /// \brief Compare UnitCellCoord
  inline bool operator>(const UnitCellCoord &A, const UnitCellCoord &B) {
    return B < A;
  }

  /// \brief Compare UnitCellCoord
  inline bool operator<=(const UnitCellCoord &A, const UnitCellCoord &B) {
    return !(A > B);
  }

  /// \brief Compare UnitCellCoord
  inline bool operator>=(const UnitCellCoord &A, const UnitCellCoord &B) {
    return !(A < B);
  }

  inline bool operator==(const UnitCellCoord &A, const UnitCellCoord &B) {
    return A.unitcell()(0) == B.unitcell()(0) &&
           A.unitcell()(1) == B.unitcell()(1) &&
           A.unitcell()(2) == B.unitcell()(2) &&
           A.sublat() == B.sublat();
  }

  /// \brief Compare UnitCellCoord
  inline bool operator!=(const UnitCellCoord &A, const UnitCellCoord &B) {
    return !(A == B);
  }


  /// \brief Print to json as [b, i, j, k]
  inline jsonParser &to_json(const UnitCellCoord &ucc_val, jsonParser &fill_json) {
    fill_json.put_array();
    fill_json.push_back(ucc_val.sublat());
    fill_json.push_back(ucc_val.unitcell()(0));
    fill_json.push_back(ucc_val.unitcell()(1));
    fill_json.push_back(ucc_val.unitcell()(2));

    return fill_json;
  }

  /// \brief Read from json [b, i, j, k]
  inline void from_json(UnitCellCoord &fill_value, const jsonParser &read_json) {

    fill_value.sublat() = read_json[0].get<Index>();
    fill_value.unitcell()(0) = read_json[1].get<Index>();
    fill_value.unitcell()(1) = read_json[2].get<Index>();
    fill_value.unitcell()(2) = read_json[3].get<Index>();

    return;
  }
}
#endif





