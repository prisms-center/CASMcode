#include <iostream>
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/SupercellEnumerator.hh"

using namespace CASM;

int main() {
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

    Lattice nigglicompare = niggli(comparelat, pg, TOL);
    Lattice nigglitest = niggli(*it, pg, TOL);

    if(nigglicompare == nigglitest) {
      std::cout << "Lattice on iteration " << l << " matched." << std::endl;
    }
    else {
      std::cout << "Lattice on iteration " << l << " did NOT match." << std::endl;
    }

    l++;
  }


  return 0;
}
