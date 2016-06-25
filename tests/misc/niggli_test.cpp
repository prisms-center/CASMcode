#include <iostream>
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/container/Counter.hh"
#include "casm/external/Eigen/Dense"
#include "casm/misc/CASM_math.hh"
#include "casm/symmetry/SymOp.hh"

using namespace CASM;

namespace Eigen {
}

namespace testing {


  void niggli_rep_test(const Lattice &in_lat) {
    int dims = 3;
    int minvol = 1;
    int maxvol = 10;

    SymGroup pg;
    in_lat.generate_point_group(pg);

    SupercellEnumerator<Lattice> latenumerator(in_lat, pg, minvol, maxvol + 1, dims);
    std::vector<Lattice> enumerated_lat(latenumerator.begin(), latenumerator.end());

    for(auto it = enumerated_lat.begin(); it != enumerated_lat.end(); ++it) {
      Lattice old_niggli = niggli(*it, pg, TOL);
      NiggliRep test_rep(old_niggli, TOL);
      NiggliRep raw_test_rep(*it, TOL);
      std::cout << "Your niggli test says " << test_rep.is_niggli() << " and " << raw_test_rep.is_niggli() << std::endl;
    }

    return;
  }

}

int main() {

  //testing::niggli(testlat);
  //testing::niggli_rep_test(Lattice::fcc());
  testing::symmetric_testing();

  return 0;
}
