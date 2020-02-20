#include "gtest/gtest.h"
#include "autotools.hh"

/// What is being tested:
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"

/// What is being used to test it:
#include <boost/filesystem/fstream.hpp>
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/SymType.hh"
#include "crystallography/TestStructures.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/app/AppIO.hh"
#include "casm/crystallography/io/VaspIO.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/crystallography/Structure.hh"

using namespace CASM;
using xtal::ScelEnumProps;
using xtal::SuperlatticeEnumerator;;

/** PRIM1 *****************************
Face-centered Cubic (FCC, cF)
1.0
0 2.0 2.0
2.0 0 2.0
2.0 2.0 0
1
D
0.00 0.00 0.00 A B C
***************************************/

void prim1_read_test(BasicStructure &struc) {

  double tol = struc.lattice().tol();

  EXPECT_EQ(almost_equal(struc.lattice()[0], Eigen::Vector3d(0.0, 2.0, 2.0), tol), true);
  EXPECT_EQ(almost_equal(struc.lattice()[1], Eigen::Vector3d(2.0, 0.0, 2.0), tol), true);
  EXPECT_EQ(almost_equal(struc.lattice()[2], Eigen::Vector3d(2.0, 2.0, 0.0), tol), true);
  EXPECT_EQ(struc.basis().size(), 1);

  // basis site 0 has three possible occupants
  EXPECT_EQ(struc.basis()[0].occupant_dof().size(), 3);

  std::string check_name[3] = {"A", "B", "C"};

  for(int i = 0; i < 3; i++) {
    // occupants are Molecule with name "A", etc.
    // Molecule are composed of AtomPosition
    // An AtomPosition 'is' a Coordinate with a Specie
    EXPECT_EQ(struc.basis()[0].occupant_dof()[i].name(), check_name[i]);
    EXPECT_EQ(almost_equal(struc.basis()[0].occupant_dof()[i].atom(0).cart(), Eigen::Vector3d(0.0, 0.0, 0.0), tol), true);
    EXPECT_EQ(struc.basis()[0].occupant_dof()[i].atom(0).name(), check_name[i]);
  }

  // FCC motif
  xtal::SymOpVector factor_group = xtal::make_factor_group(struc);
  EXPECT_EQ(48, factor_group.size());
}

/** PRIM2 *****************************
Face-centered Cubic (FCC, cF)
1.0
4.0 0.0 0.0
0.0 4.0 0.0
0.0 0.0 4.0
2 1 1
D
0.00 0.00 0.00 A B C :: A
0.5 0.5 0.0 A B C :: A
0.5 0.0 0.5 A B C :: B
0.0 0.5 0.5 A B C :: C
***************************************/

void prim2_read_test(BasicStructure &struc) {

  double tol = struc.lattice().tol();

  EXPECT_EQ(almost_equal(struc.lattice()[0], Eigen::Vector3d(4.0, 0.0, 0.0), tol), true);
  EXPECT_EQ(almost_equal(struc.lattice()[1], Eigen::Vector3d(0.0, 4.0, 0.0), tol), true);
  EXPECT_EQ(almost_equal(struc.lattice()[2], Eigen::Vector3d(0.0, 0.0, 4.0), tol), true);
  EXPECT_EQ(struc.basis().size(), 4);

  // basis site 0 has three possible occupants
  EXPECT_EQ(struc.basis()[0].occupant_dof().size(), 3);

  std::string check_name[3] = {"A", "B", "C"};
  int check_value[4] = {0, 0, 1, 2};

  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 3; j++) {
      // occupants are Molecule with name "A", etc.
      // Molecule are composed of AtomPosition
      // An AtomPosition 'is' a Coordinate with a Specie
      EXPECT_EQ(struc.basis()[i].occupant_dof()[j].name(), check_name[j]);
      EXPECT_EQ(almost_equal(struc.basis()[i].occupant_dof()[j].atom(0).cart(), Eigen::Vector3d(0.0, 0.0, 0.0), tol), true);
      EXPECT_EQ(struc.basis()[i].occupant_dof()[j].atom(0).name(), check_name[j]);
    }
    /* EXPECT_EQ(struc.basis()[i].occupant_dof().value(), check_value[i]); */
  }


  //Modify the structure that there's different occupants at each site
  std::vector<Site> new_basis;
  new_basis.emplace_back(struc.basis()[0], "A");
  new_basis.emplace_back(struc.basis()[1], "A");
  new_basis.emplace_back(struc.basis()[2], "B");
  new_basis.emplace_back(struc.basis()[3], "C");
  struc.set_basis(new_basis);


  // ordering on FCC motif
  xtal::SymOpVector factor_group = xtal::make_factor_group(struc);
  EXPECT_EQ(16, factor_group.size());
}


/** POS1 *****************************
Face-centered Cubic (FCC, cF)
1.0
4.0 0.0 0.0
0.0 4.0 0.0
0.0 0.0 4.0
A B C
2 1 1
D
0.00 0.00 0.00 A
0.5 0.5 0.0 A
0.5 0.0 0.5 B
0.0 0.5 0.5 C
***************************************/


/** POS1_vasp5_out ********************
Face-centered Cubic (FCC, cF)
 1.00000000
       4.00000000      0.00000000      0.00000000
       0.00000000      4.00000000      0.00000000
       0.00000000      0.00000000      4.00000000
A B C
2 1 1
Direct
   0.0000000   0.0000000   0.0000000
   0.5000000   0.5000000   0.0000000
   0.5000000   0.0000000   0.5000000
   0.0000000   0.5000000   0.5000000
***************************************/


void pos1_read_test(BasicStructure &struc) {

  double tol = struc.lattice().tol();

  EXPECT_EQ(almost_equal(struc.lattice()[0], Eigen::Vector3d(4.0, 0.0, 0.0), tol), true);
  EXPECT_EQ(almost_equal(struc.lattice()[1], Eigen::Vector3d(0.0, 4.0, 0.0), tol), true);
  EXPECT_EQ(almost_equal(struc.lattice()[2], Eigen::Vector3d(0.0, 0.0, 4.0), tol), true);
  EXPECT_EQ(struc.basis().size(), 4);

  std::string check_name[4] = {"A", "A", "B", "C"};

  for(int i = 0; i < 4; i++) {
    // basis site 0 and 1 have one possible occupant
    EXPECT_EQ(struc.basis()[i].occupant_dof().size(), 1);

    // occupants are Molecule with name "A", etc.
    // Molecule are composed of AtomPosition
    // An AtomPosition 'is' a Coordinate with a Specie
    EXPECT_EQ(struc.basis()[i].occupant_dof()[0].name(), check_name[i]);
    EXPECT_EQ(almost_equal(struc.basis()[i].occupant_dof()[0].atom(0).cart(), Eigen::Vector3d(0.0, 0.0, 0.0), tol), true);
    EXPECT_EQ(struc.basis()[i].occupant_dof()[0].atom(0).name(), check_name[i]);
  }

  // FCC structure
  xtal::SymOpVector factor_group = xtal::make_factor_group(struc);
  EXPECT_EQ(16, factor_group.size());
}

namespace {
  fs::path crystallography_test_directory() {
    return autotools::abs_srcdir() + "/tests/unit/crystallography";
  }
}


TEST(BasicStructureSiteTest, PRIM1Test) {

  fs::path testdir = ::crystallography_test_directory();

  // Read in test PRIM and run tests
  BasicStructure struc(fs::path(testdir / "PRIM1.txt"));
  prim1_read_test(struc);

  // Write test PRIM back out
  fs::path tmp_file = testdir / "PRIM1_out.txt";
  write_prim(struc, tmp_file, FRAC);

  // Read new file and run tests again
  BasicStructure struc2(read_prim(tmp_file, TOL));
  prim1_read_test(struc2);

}

TEST(BasicStructureSiteTest, PRIM2Test) {

  fs::path testdir = ::crystallography_test_directory();

  // Read in test PRIM and run tests
  BasicStructure struc(fs::path(testdir / "PRIM2.txt"));
  prim2_read_test(struc);
}

TEST(BasicStructureSiteTest, PRIM3Test) {

  fs::path testdir = ::crystallography_test_directory();

  // Read in an incorrectly formatted PRIM and check that an exception is thrown
  EXPECT_THROW(BasicStructure(fs::path(testdir / "PRIM3.txt")), std::runtime_error);

}

TEST(BasicStructureSiteTest, POS1Test) {

  fs::path testdir = ::crystallography_test_directory();

  // Read in test PRIM and run tests
  BasicStructure struc(fs::path(testdir / "POS1.txt"));
  pos1_read_test(struc);

  // Write test PRIM back out
  fs::path tmp_file = testdir / "POS1_out.txt";
  fs::ofstream sout(tmp_file);
  VaspIO::PrintPOSCAR printer(make_simple_structure(struc), struc.title());
  printer.set_append_atom_names_off();
  printer.print(sout);
  sout.close();

  // Read new file and run tests again
  BasicStructure struc2(fs::path(testdir / "POS1_out.txt"));
  pos1_read_test(struc2);

}

TEST(BasicStructureSiteTest, POS1Vasp5Test) {

  fs::path testdir = ::crystallography_test_directory();

  // Read in test PRIM and run tests
  BasicStructure struc(fs::path(testdir / "POS1.txt"));
  pos1_read_test(struc);

  // Write test PRIM back out
  fs::path tmp_file = testdir / "POS1_vasp5_out.txt";
  fs::ofstream sout(tmp_file);
  VaspIO::PrintPOSCAR(make_simple_structure(struc), struc.title()).print(sout);
  sout.close();

  // Read new file and run tests again
  BasicStructure struc2(fs::path(testdir / "POS1_vasp5_out.txt"));
  pos1_read_test(struc2);

}

TEST(BasicStructureSiteTest, IsPrimitiveTest) {

  Structure prim(test::ZrO_prim());

  const SymGroup effective_pg = prim.factor_group();
  std::cout << effective_pg.size() << std::endl;

  ScelEnumProps enum_props(1, 7);
  SuperlatticeEnumerator scel_enum(effective_pg.begin(), effective_pg.end(), prim.lattice(), enum_props);

  for(auto it = scel_enum.begin(); it != scel_enum.end(); ++it) {

    Structure super = prim.create_superstruc(*it);
    EXPECT_EQ(super.lattice().is_right_handed(), true);

    Structure new_prim = Structure(xtal::make_primitive(super));

    auto is_trans_pair = xtal::is_superlattice(super.lattice(), prim.lattice(), TOL);

    EXPECT_EQ(new_prim.lattice().is_right_handed(), true);
    EXPECT_EQ(new_prim.lattice().is_right_handed(), super.lattice().is_right_handed());
  }

}
