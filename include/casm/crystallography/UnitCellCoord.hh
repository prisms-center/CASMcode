#ifndef UNITCELLCOORD_HH
#define UNITCELLCOORD_HH

#include <iostream>

namespace CASM {
  /**
   * This class is a lightweight version of Site. Given 4
   * indices b, i, j, k it is possible to infer what in
   * a structure you're referring to.
   *
   * b represents the basis atom in the prim
   * i, j, k tell you which primitive volume in the
   * superstructure you're in.
   */

  class UnitCellCoord {
  public:
    long int coord[4];	// b i j k since the supercell maintains this ordering

    ///Default construction setting all values to zero
    UnitCellCoord();

    ///Explicit construction in order b, i, j, k
    UnitCellCoord(long int i0, long int i1, long int i2, long int i3);

    long int &operator[](int i1);

    const long int &operator[](int i1) const;

    bool operator==(UnitCellCoord B) const;

    bool operator!=(UnitCellCoord B) const;

    bool operator<(UnitCellCoord B) const;

    UnitCellCoord operator+(UnitCellCoord B) const;

    UnitCellCoord operator-(UnitCellCoord B) const;        //John G changed for bijk

    UnitCellCoord operator*(long int i) const;         //John G changed for bijk

    ///Allow << stream operator for printing
    friend std::ostream &operator<< (std::ostream &stream, const UnitCellCoord &tucc);

    ///Explicitly redefine all values in order b, i, j, k
    void set(long int i0, long int i1, long int i2, long int i3);
  };

  class jsonParser;

  jsonParser &to_json(const UnitCellCoord &ucc_val, jsonParser &fill_json);
  void from_json(UnitCellCoord &fill_value, const jsonParser &read_json);

}
#endif


