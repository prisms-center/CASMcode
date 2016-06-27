#include <iostream>
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/container/Counter.hh"
#include "casm/external/Eigen/Dense"
#include "casm/misc/CASM_math.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/casm_io/VaspIO.hh"
#include "casm/clex/Supercell.hh"


using namespace CASM;

namespace Eigen {
}

namespace testing {


  /*
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
  */

  void obvious_check() {
    Structure non_reduced("./not_reduced.vasp");
    Lattice non_reduced_lat = non_reduced.lattice();

    Structure reduced("./reduced.vasp");
    Lattice reduced_lat = reduced.lattice();

    std::cout << "The non reduced structure appears to be niggli? " << is_niggli(non_reduced_lat, CASM::TOL) << std::endl;
    std::cout << "The reduced structure appears to be niggli? " << is_niggli(reduced_lat, CASM::TOL) << std::endl;

    Lattice made_niggli = niggli(reduced_lat, CASM::TOL);
    Structure made_niggli_struc = reduced.create_superstruc(made_niggli);

    std::cout << "DEBUGGING: made_niggli.lat_column_mat().determinant() is " << made_niggli.lat_column_mat().determinant() << std::endl;


    VaspIO::PrintPOSCAR printer(reduced_lat);
    printer.set_atom_names_on();
    printer.set_coord_mode(CASM::FRAC);
    printer.print(std::cout);

    NiggliRep nigrep(made_niggli);
    std::cout << nigrep.meets_criteria_1(CASM::TOL) << std::endl;
    std::cout << nigrep.meets_criteria_2(CASM::TOL) << std::endl;
    std::cout << nigrep.meets_criteria_3(CASM::TOL) << std::endl;
    std::cout << nigrep.meets_criteria_4(CASM::TOL) << std::endl;
    std::cout << nigrep.meets_criteria_5(CASM::TOL) << std::endl;
    std::cout << nigrep.meets_criteria_6(CASM::TOL) << std::endl;
    std::cout << nigrep.meets_criteria_7(CASM::TOL) << std::endl;
    std::cout << nigrep.meets_criteria_8(CASM::TOL) << std::endl;
    std::cout << "DEBUGGING: nigrep.A() is " << nigrep.A() << std::endl;
    std::cout << "DEBUGGING: nigrep.B() is " << nigrep.B() << std::endl;
    std::cout << "DEBUGGING: nigrep.C() is " << nigrep.C() << std::endl;
    std::cout << "DEBUGGING: nigrep.ksi() is " << nigrep.ksi() << std::endl;
    std::cout << "DEBUGGING: nigrep.eta() is " << nigrep.eta() << std::endl;
    std::cout << "DEBUGGING: nigrep.zeta() is " << nigrep.zeta() << std::endl;
    double clobber = nigrep.A() + nigrep.B() + nigrep.C() + nigrep.ksi() + nigrep.eta() + nigrep.zeta();
    std::cout << "DEBUGGING: clobber is " << clobber << std::endl;



    std::cout << std::endl;

    std::cout << nigrep.metrical_matrix() << std::endl;
    std::cout << std::endl;
    std::cout << reduced_lat.lat_column_mat().inverse()*made_niggli.lat_column_mat() << std::endl;

    return;
  }

  void old_enumerations_test(const std::string &namelist) {
    std::ifstream path_stream(namelist.c_str());

    std::string path_str;
    std::string prim_path_str;

    std::getline(path_stream, prim_path_str);

    SymGroup old_primitive_point;
    Structure old_primitive_struc(prim_path_str);
    old_primitive_struc.lattice().generate_point_group(old_primitive_point);

    SymGroup new_primitive_point;
    Structure new_primitive_struc(prim_path_str);
    new_primitive_struc.lattice().generate_point_group(new_primitive_point);

    path_stream.close();
    path_stream.open(namelist.c_str());

    while(std::getline(path_stream, path_str)) {

      Structure old_struc(path_str);
      NiggliRep old_rep(old_struc.lattice());

      if(!old_rep.is_niggli_type2(CASM::TOL) && !old_rep.is_niggli_type1(CASM::TOL)) {
        std::cout << path_str << " was never Niggli to begin with!!" << std::endl;
      }

      Lattice old_canonical_lat = old_struc.lattice();
      Lattice new_canonical_lat = canonical_equivalent_lattice(old_canonical_lat, new_primitive_point, CASM::TOL);

      Eigen::Matrix3d old_trans_mat, new_trans_mat;

      assert(old_canonical_lat.is_supercell_of(old_primitive_struc.lattice(), old_trans_mat));
      assert(new_canonical_lat.is_supercell_of(new_primitive_struc.lattice(), new_trans_mat));

      std::string old_name_str, new_name_str;
      old_name_str = generate_name(old_trans_mat.cast<int>());
      new_name_str = generate_name(new_trans_mat.cast<int>());

      if(old_name_str != new_name_str) {
        std::cout << "Name has changed: " << old_name_str << " -> " << new_name_str << std::endl;
      }

      if(!(new_canonical_lat == old_canonical_lat)) {
        std::cout << path_str << " is different now..." << std::endl;

        Eigen::Matrix3d old_to_new_trans = Eigen::Matrix3d::Zero();
        bool transformable;

        transformable = new_canonical_lat.is_supercell_of(old_canonical_lat, new_primitive_point, old_to_new_trans);

        if(transformable) {
          std::cout << "Can go from old to new with:" << std::endl;
          std::cout << old_to_new_trans.cast<int>() << std::endl << std::endl;
        }

        else {
          std::cout << "Could not transform" << std::endl;
        }
      }


    }
    return;
  }

}

int main() {

  //testing::niggli(testlat);
  //testing::niggli_rep_test(Lattice::fcc());
  testing::obvious_check();
  testing::old_enumerations_test("enumerated.txt");

  return 0;
}
