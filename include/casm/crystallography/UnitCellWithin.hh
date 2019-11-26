#ifndef UNITCELLWITHIN_HH
#define UNITCELLWITHIN_HH

#include "casm/external/Eigen/Core"

namespace CASM {
  namespace xtal {
    class UnitCell;

    /**
     * Handles bringing a UnitCell (i,j,k values) within a particular
     * superlattice. Provide the tiling unit and superlattice at construction,
     * and use this class to bring any external UnitCell within the
     * superlattice
     */

    struct UnitCellWithin {
      UnitCell operator()(const UnitCell &ijk) const;

    private:
      typedef Eigen::Matrix<long, 3, 3> matrix_type;
      typedef Eigen::Matrix<long, 3, 1> vector_type;
    };
  }
}

#endif
