#ifndef SYMBASISPERMUTE_HH
#define SYMBASISPERMUTE_HH

#include <iostream>
#include <cmath>

namespace CASM {
  class MasterSymGroup;
  class Lattice;


  ///\brief  SymBasisPermute describes how a symmetry operation permutes atoms in a basis
  /// It is represented as an Array<UnitCellCoord> that describes how
  /// site (b,0,0,0) -> (b_new,i,j,k) under application of symmetry
  class SymBasisPermute: public SymOpRepresentation {
  private:

    /// Array of UnitCellCoords, of length basis.size(). Site (b,0,0,0) goes to
    /// (b_new,i,j,k)==ucc_permute[b] after application of symmetry.
    Array< UnitCellCoord > ucc_permute;


  public:

    /// Initialize a SymPermutation with the permutation array.
    SymBasisPermute(const Array<UnitCellCoord> &init_permute) : ucc_permute(init_permute) {

    };

    /// Return pointer to a copy of this SymBasisPermute
    SymOpRepresentation *copy() const {
      return new SymBasisPermute(*this);
    };

    /// Access the permutation array 'ucc_permute'
    Array<UnitCellCoord> const *get_ucc_permutation() const {
      return &ucc_permute;
    };

    jsonParser &to_json(jsonParser &json) const {
      return json;
    };
    void from_json(const jsonParser &json) {};

  };


}
#endif
