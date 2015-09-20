#include "casm/crystallography/UnitCellCoord.hh"

#include <iostream>
#include "casm/casm_io/jsonParser.hh"

namespace CASM {
  UnitCellCoord::UnitCellCoord() {
    coord[0] = 0;
    coord[1] = 0;
    coord[2] = 0;
    coord[3] = 0;

  }

  UnitCellCoord::UnitCellCoord(long int i0, long int i1, long int i2, long int i3) {
    coord[0] = i0;
    coord[1] = i1;
    coord[2] = i2;
    coord[3] = i3;

  }

  long int &UnitCellCoord::operator[](int i1) {
    return coord[i1];
  }

  const long int &UnitCellCoord::operator[](int i1) const {
    return coord[i1];
  }

  bool UnitCellCoord::operator==(UnitCellCoord B) const {
    if(coord[0] == B[0])
      if(coord[1] == B[1])
        if(coord[2] == B[2])
          if(coord[3] == B[3]) {
            return 1;
          }

    return 0;
  }

  bool UnitCellCoord::operator!=(UnitCellCoord B) const {
    if(coord[0] == B[0])
      if(coord[1] == B[1])
        if(coord[2] == B[2])
          if(coord[3] == B[3]) {
            return 0;
          }

    return 1;
  }

  bool UnitCellCoord::operator<(UnitCellCoord B) const {       //John G changed for bijk
    if(coord[1] < B[1]) {
      return true;
    }
    else if(coord[2] == B[2]) {
      if(coord[2] < B[2]) {
        return true;
      }
      else if(coord[2] == B[2]) {
        if(coord[3] < B[3]) {
          return true;
        }

      }
    }

    return false;
  }

  //UnitCellCoord operator+(UnitCellCoord B) {
  //  UnitCellCoord A(*this);

  //  A[0] =  B[0];
  //  A[1] += B[1];
  //  A[2] += B[2];
  //  A[3] += B[3];

  //  return A;
  //}

  UnitCellCoord UnitCellCoord::operator+(UnitCellCoord B) const {
    UnitCellCoord A(*this);

    A[0] =  B[0]; //Is this really what we want? Adding UCC for different basis seems weird
    A[1] += B[1];
    A[2] += B[2];
    A[3] += B[3];

    return A;
  }

  //UnitCellCoord operator-(UnitCellCoord B) {       //John G changed for bijk
  //  UnitCellCoord A(*this);

  //  A[1] -= B[1];
  //  A[2] -= B[2];
  //  A[3] -= B[3];
  //  //A[0] unchanged

  //  return A;
  //}

  UnitCellCoord UnitCellCoord::operator-(UnitCellCoord B) const {       //John G changed for bijk
    UnitCellCoord A(*this);

    A[1] -= B[1];
    A[2] -= B[2];
    A[3] -= B[3];
    //A[0] unchanged

    return A;
  }

  UnitCellCoord UnitCellCoord::operator*(long int i) const {       //John G changed for bijk
    UnitCellCoord A(*this);

    A[1] *= i;
    A[2] *= i;
    A[3] *= i;

    return A;
  }

  std::ostream &operator<<(std::ostream &stream, const UnitCellCoord &tucc) {
    stream << tucc[0] << " " << tucc[1] << " " << tucc[2] << " " << tucc[3] << "\n";
    return stream;
  }

  // friend istream &operator>>(istream &instream, UnitCellCoord &c) {
  //   instream >> c[0] >> c[1] >> c[2] >> c[3] ;
  //   return instream;
  // }

  void UnitCellCoord::set(long int i0, long int i1, long int i2, long int i3) {
    coord[0] = i0;
    coord[1] = i1;
    coord[2] = i2;
    coord[3] = i3;
    return;
  }


  //Can't put this in UntiCellCoord file unless yous split it into hh and cc, which apparently we decided against
  jsonParser &to_json(const UnitCellCoord &ucc_val, jsonParser &fill_json) {
    fill_json.put_array();
    fill_json.push_back(ucc_val[0]);
    fill_json.push_back(ucc_val[1]);
    fill_json.push_back(ucc_val[2]);
    fill_json.push_back(ucc_val[3]);

    return fill_json;
  }

  void from_json(UnitCellCoord &fill_value, const jsonParser &read_json) {
    long int bb, ii, jj, kk;
    bb = read_json[0].get<long int>();
    ii = read_json[1].get<long int>();
    jj = read_json[2].get<long int>();
    kk = read_json[3].get<long int>();

    fill_value.set(bb, ii, jj, kk);
    //fill_value.set(read_json["b"].get<long int>(), read_json["i"].get<long int>(), read_json["j"].get<long int>(), read_json["k"].get<long int>());
    return;
  }


}

