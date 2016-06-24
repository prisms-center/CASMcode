#include <iostream>
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/container/Counter.hh"
#include "casm/external/Eigen/Dense"
#include "casm/misc/CASM_math.hh"
#include "casm/symmetry/SymOp.hh"

using namespace CASM;

namespace Eigen {
}

namespace testing {

  void single_dimension_test() {
    Lattice testlat = Lattice::fcc();
    SymGroup pg;
    testlat.generate_point_group(pg);

    int dims = 1;
    int minvol = 1;
    int maxvol = 10;

    SupercellEnumerator<Lattice> latenumerator(testlat, pg, minvol, maxvol + 1, dims);
    std::vector<Lattice> enumerated_lat(latenumerator.begin(), latenumerator.end());

    std::cout << "Enumerated from " << minvol << " to " << maxvol << " and got " << enumerated_lat.size() << " lattices" << std::endl;


    int l = 1;
    for(auto it = enumerated_lat.begin(); it != enumerated_lat.end(); ++it) {
      Eigen::Matrix3i comp_transmat;
      comp_transmat << (l), 0, 0,
                    0, 1, 0,
                    0, 0, 1;

      Lattice comparelat = make_supercell(testlat, comp_transmat);

      //Lattice nigglicompare = niggli(comparelat, TOL);
      //Lattice nigglitest = niggli(*it, TOL);

      //Lattice nigglicompare = niggli(comparelat, pg, TOL);
      //Lattice nigglitest = niggli(*it, pg, TOL);


      Lattice nigglicompare = canonical_equivalent_lattice(comparelat, pg, TOL);
      Lattice nigglitest = canonical_equivalent_lattice(*it, pg, TOL);


      if(nigglicompare == nigglitest)
        //if(nigglicompare.is_equivalent(nigglitest, CASM::TOL))
      {
        std::cout << "Lattice on iteration " << l << " matched." << std::endl;
      }
      else {
        std::cout << "Lattice on iteration " << l << " did NOT match." << std::endl;
      }

      l++;
    }
    return;
  }

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

  void symmetric_testing() {
    Eigen::MatrixXd symmat(5, 5), persymmat(4, 4), nonsymmat;

    symmat << 1, 2, 3, 4, 5,
           2, 6, 7, 8, 9,
           3, 7, 10, 11, 12,
           4, 8, 11, 13, 14,
           5, 9, 12, 14, 15;

    std::cout << "symmetric matrix is symmetric? " << is_symmetric(symmat) << std::endl;
    std::cout << "symmetric matrix is persymmetric? " << is_persymmetric(symmat) << std::endl;

    persymmat << 4, 3, 2, 1,
              7, 6, 5, 2,
              9, 8, 6, 3,
              10, 9, 7, 4;

    std::cout << "persymmetric matrix is symmetric? " << is_symmetric(persymmat) << std::endl;
    std::cout << "persymmetric matrix is persymmetric? " << is_persymmetric(persymmat) << std::endl;

    return;
  }

}

int main() {

  //testing::niggli(testlat);
  testing::single_dimension_test();
  //testing::niggli_rep_test(Lattice::fcc());
  testing::symmetric_testing();

  return 0;
}
