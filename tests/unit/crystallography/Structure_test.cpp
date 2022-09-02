#include "Common.hh"
#include "gtest/gtest.h"

/// What is being tested:
#include "casm/crystallography/Structure.hh"

/// What is being used to test it:
#include <boost/filesystem/fstream.hpp>

#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/crystallography/io/VaspIO.hh"
#include "casm/misc/CASM_Eigen_math.hh"

using namespace CASM;

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

void prim1_read_test(Structure &struc) {
  double tol = 1e-5;

  EXPECT_EQ(
      almost_equal(struc.lattice()[0], Eigen::Vector3d(0.0, 2.0, 2.0), tol),
      true);
  EXPECT_EQ(
      almost_equal(struc.lattice()[1], Eigen::Vector3d(2.0, 0.0, 2.0), tol),
      true);
  EXPECT_EQ(
      almost_equal(struc.lattice()[2], Eigen::Vector3d(2.0, 2.0, 0.0), tol),
      true);
  EXPECT_EQ(struc.basis().size(), 1);

  // basis site 0 has three possible occupants
  EXPECT_EQ(struc.basis()[0].occupant_dof().size(), 3);

  std::string check_name[3] = {"A", "B", "C"};

  for (int i = 0; i < 3; i++) {
    // occupants are Molecule with name "A", etc.
    // Molecule are composed of AtomPosition
    // An AtomPosition 'is' a Coordinate with a Specie
    EXPECT_EQ(struc.basis()[0].occupant_dof()[i].name(), check_name[i]);
    EXPECT_EQ(almost_equal(struc.basis()[0].occupant_dof()[i].atom(0).cart(),
                           Eigen::Vector3d(0.0, 0.0, 0.0), tol),
              true);
    EXPECT_EQ(struc.basis()[0].occupant_dof()[i].atom(0).name(), check_name[i]);
  }

  // FCC motif
  EXPECT_EQ(48, struc.factor_group().size());
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

void prim2_read_test(Structure &struc) {
  double tol = 1e-5;

  EXPECT_EQ(
      almost_equal(struc.lattice()[0], Eigen::Vector3d(4.0, 0.0, 0.0), tol),
      true);
  EXPECT_EQ(
      almost_equal(struc.lattice()[1], Eigen::Vector3d(0.0, 4.0, 0.0), tol),
      true);
  EXPECT_EQ(
      almost_equal(struc.lattice()[2], Eigen::Vector3d(0.0, 0.0, 4.0), tol),
      true);
  EXPECT_EQ(struc.basis().size(), 4);

  // basis site 0 has three possible occupants
  EXPECT_EQ(struc.basis()[0].occupant_dof().size(), 3);

  std::string check_name[3] = {"A", "B", "C"};
  int check_value[4] = {0, 0, 1, 2};

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++) {
      // occupants are Molecule with name "A", etc.
      // Molecule are composed of AtomPosition
      // An AtomPosition 'is' a Coordinate with a Specie
      EXPECT_EQ(struc.basis()[i].occupant_dof()[j].name(), check_name[j]);
      EXPECT_EQ(almost_equal(struc.basis()[0].occupant_dof()[j].atom(0).cart(),
                             Eigen::Vector3d(0.0, 0.0, 0.0), tol),
                true);
      EXPECT_EQ(struc.basis()[i].occupant_dof()[j].atom(0).name(),
                check_name[j]);
    }
  }

  // Modify the structure that there's different occupants at each site
  std::vector<xtal::Site> new_basis;
  new_basis.emplace_back(struc.basis()[0], "A");
  new_basis.emplace_back(struc.basis()[1], "A");
  new_basis.emplace_back(struc.basis()[2], "B");
  new_basis.emplace_back(struc.basis()[3], "C");

  xtal::BasicStructure struc_with_new_basis = struc.structure();
  struc_with_new_basis.set_basis(new_basis);
  struc = Structure(struc_with_new_basis);

  // ordering on FCC motif
  EXPECT_EQ(16, struc.factor_group().size());
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

void pos1_read_test(Structure &struc) {
  double tol = 1e-5;

  EXPECT_EQ(
      almost_equal(struc.lattice()[0], Eigen::Vector3d(4.0, 0.0, 0.0), tol),
      true);
  EXPECT_EQ(
      almost_equal(struc.lattice()[1], Eigen::Vector3d(0.0, 4.0, 0.0), tol),
      true);
  EXPECT_EQ(
      almost_equal(struc.lattice()[2], Eigen::Vector3d(0.0, 0.0, 4.0), tol),
      true);
  EXPECT_EQ(struc.basis().size(), 4);

  std::string check_name[4] = {"A", "A", "B", "C"};

  for (int i = 0; i < 4; i++) {
    // basis site 0 and 1 have one possible occupant
    EXPECT_EQ(struc.basis()[i].occupant_dof().size(), 1);

    // occupants are Molecule with name "A", etc.
    // Molecule are composed of AtomPosition
    // An AtomPosition 'is' a Coordinate with a Specie
    EXPECT_EQ(struc.basis()[i].occupant_dof()[0].name(), check_name[i]);
    EXPECT_EQ(almost_equal(struc.basis()[i].occupant_dof()[0].atom(0).cart(),
                           Eigen::Vector3d(0.0, 0.0, 0.0), tol),
              true);
    EXPECT_EQ(struc.basis()[i].occupant_dof()[0].atom(0).name(), check_name[i]);
  }

  // FCC structure
  EXPECT_EQ(16, struc.factor_group().size());
}

class StructureTest : public testing::Test {
 protected:
  fs::path datadir;
  test::TmpDir tmpdir;

  StructureTest() : datadir(test::data_dir("crystallography")) {}

  // // Can use this to check write failures:
  // void TearDown() {
  //   if(HasFailure()) {
  //     std::cout << "tmpdir: " << tmpdir.path() << std::endl;
  //     tmpdir.do_not_remove_on_destruction();
  //   }
  // }
};

TEST_F(StructureTest, PRIM1Test) {
  // Read in test PRIM and run tests
  Structure struc(datadir / "PRIM1.txt");
  prim1_read_test(struc);

  // Write test PRIM back out
  fs::path tmp_file = tmpdir.path() / "PRIM1_out.txt";
  write_prim(struc, tmp_file, FRAC);

  // Read new file and run tests again
  Structure struc2(read_prim(tmp_file, TOL));
  prim1_read_test(struc2);
}

TEST_F(StructureTest, PRIM2Test) {
  // Read in test PRIM and run tests
  Structure struc(datadir / "PRIM2.txt");
  prim2_read_test(struc);
}

TEST_F(StructureTest, POS1Test) {
  // Read in test PRIM and run tests
  Structure struc(datadir / "POS1.txt");
  pos1_read_test(struc);

  // Write test PRIM back out
  fs::path tmp_file = tmpdir.path() / "POS1_out.txt";
  fs::ofstream sout(tmp_file);
  VaspIO::PrintPOSCAR printer(xtal::make_simple_structure(struc),
                              struc.structure().title());
  printer.set_append_atom_names_off();
  printer.print(sout);
  sout.close();

  // Read new file and run tests again
  Structure struc2(tmp_file);
  pos1_read_test(struc2);
}

TEST_F(StructureTest, POS1Vasp5Test) {
  // Read in test PRIM and run tests
  Structure struc(datadir / "POS1.txt");
  pos1_read_test(struc);

  // Write test PRIM back out
  fs::path tmp_file = tmpdir.path() / "POS1_vasp5_out.txt";
  fs::ofstream sout(tmp_file);
  VaspIO::PrintPOSCAR(xtal::make_simple_structure(struc),
                      struc.structure().title())
      .print(sout);
  sout.close();

  // Read new file and run tests again
  Structure struc2(tmp_file);
  pos1_read_test(struc2);
}

TEST_F(StructureTest, POS1jsonPrimTest) {
  // Read in test PRIM and run tests
  Structure struc(datadir / "POS1.txt");
  pos1_read_test(struc);

  // Write test PRIM back out
  fs::path tmp_file = tmpdir.path() / "POS1_prim.json";
  jsonParser json;
  write_prim(struc, json, FRAC);
  fs::ofstream sout(tmp_file);
  json.print(sout);
  sout.close();

  // Read new file and run tests again
  struc = Structure(read_prim(tmp_file, TOL));
  pos1_read_test(struc);
}
