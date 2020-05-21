#include "autotools.hh"
#include "gtest/gtest.h"

/// What is being tested:
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"

/// What is being used to test it:
#include "casm/app/AppIO.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/crystallography/io/VaspIO.hh"
#include "casm/external/Eigen/src/Core/Matrix.h"
#include "casm/misc/CASM_Eigen_math.hh"
#include "crystallography/TestStructures.hh"
#include <boost/filesystem/fstream.hpp>
#include <fstream>

using namespace CASM;
using xtal::ScelEnumProps;
using xtal::SuperlatticeEnumerator;
;

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

  // Modify the structure that there's different occupants at each site
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

  BasicStructure read_structure(const fs::path &poscar_path) {
    std::ifstream poscar_stream(poscar_path.string());
    return BasicStructure::from_poscar_stream(poscar_stream);
  }
} // namespace

TEST(BasicStructureSiteTest, PRIM1Test) {
  fs::path testdir = ::crystallography_test_directory();

  // Read in test PRIM and run tests
  BasicStructure struc = ::read_structure(fs::path(testdir / "PRIM1.txt"));
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
  BasicStructure struc = ::read_structure(fs::path(testdir / "PRIM2.txt"));
  prim2_read_test(struc);
}

TEST(BasicStructureSiteTest, PRIM3Test) {

  fs::path testdir = ::crystallography_test_directory();

  // Read in an incorrectly formatted PRIM and check that an exception is thrown
  EXPECT_THROW(::read_structure(fs::path(testdir / "PRIM3.txt")), std::runtime_error);
}

TEST(BasicStructureSiteTest, POS1Test) {

  fs::path testdir = ::crystallography_test_directory();

  // Read in test PRIM and run tests
  BasicStructure struc = ::read_structure(fs::path(testdir / "POS1.txt"));
  pos1_read_test(struc);

  // Write test PRIM back out
  fs::path tmp_file = testdir / "POS1_out.txt";
  fs::ofstream sout(tmp_file);
  VaspIO::PrintPOSCAR printer(make_simple_structure(struc), struc.title());
  printer.set_append_atom_names_off();
  printer.print(sout);
  sout.close();

  // Read new file and run tests again
  BasicStructure struc2 = ::read_structure(fs::path(testdir / "POS1_out.txt"));
  pos1_read_test(struc2);
}

TEST(BasicStructureSiteTest, POS1Vasp5Test) {

  fs::path testdir = ::crystallography_test_directory();

  // Read in test PRIM and run tests
  BasicStructure struc = ::read_structure(fs::path(testdir / "POS1.txt"));
  pos1_read_test(struc);

  // Write test PRIM back out
  fs::path tmp_file = testdir / "POS1_vasp5_out.txt";
  fs::ofstream sout(tmp_file);
  VaspIO::PrintPOSCAR(make_simple_structure(struc), struc.title()).print(sout);
  sout.close();

  // Read new file and run tests again
  BasicStructure struc2 = ::read_structure(fs::path(testdir / "POS1_vasp5_out.txt"));
  pos1_read_test(struc2);
}

TEST(BasicStructureSiteTest, MakeSuperstructure) {
  Eigen::Matrix3i transf_mat;
  transf_mat << 2, 0, 0, 0, 2, 0, 0, 0, 2;

  BasicStructure prim = test::ZrO_prim();
  BasicStructure superstruc = xtal::make_superstructure(prim, transf_mat);

  EXPECT_EQ(superstruc.basis().size(), transf_mat.determinant() * prim.basis().size());

  for(const Site &prim_site : prim.basis()) {
    int match_count = 0;
    for(const Site &super_site : superstruc.basis()) {
      if(super_site.compare(prim_site)) {
        ++match_count;
      }
    }
    EXPECT_EQ(match_count, 1);
  }
}

TEST(BasicStructureSiteTest, IsPrimitiveTest) {

  Structure prim(test::ZrO_prim());

  const SymGroup effective_pg = prim.factor_group();
  std::cout << effective_pg.size() << std::endl;

  ScelEnumProps enum_props(1, 7);
  SuperlatticeEnumerator scel_enum(effective_pg.begin(), effective_pg.end(), prim.lattice(), enum_props);

  for(auto it = scel_enum.begin(); it != scel_enum.end(); ++it) {
    Eigen::Matrix3l transformation_matix = xtal::make_transformation_matrix_to_super(prim.lattice(), *it, prim.lattice().tol());
    BasicStructure super = xtal::make_superstructure(prim.lattice(), transformation_matix);

    EXPECT_EQ(super.lattice().is_right_handed(), true);

    Structure new_prim = Structure(xtal::make_primitive(super));

    auto is_trans_pair = xtal::is_superlattice(super.lattice(), prim.lattice(), TOL);

    EXPECT_EQ(new_prim.lattice().is_right_handed(), true);
    EXPECT_EQ(new_prim.lattice().is_right_handed(), super.lattice().is_right_handed());
  }
}

//************************************************************************//

class HexagonalSuperStructureTest : public testing::Test {
protected:
  void SetUp() override {
    fs::path hcp_stack_file = ::crystallography_test_directory() / "hcp_stack3.vasp";
    std::ifstream hcp_stack_stream(hcp_stack_file.string());
    xtal::BasicStructure hcp_stack_structure = xtal::BasicStructure::from_poscar_stream(hcp_stack_stream);

    hcp_3stack_ptr.reset(new xtal::BasicStructure(std::move(hcp_stack_structure)));
    hcp_primitive_ptr.reset(new xtal::BasicStructure(xtal::make_primitive(*hcp_3stack_ptr)));


    fs::path hcp_read_prim_path = ::crystallography_test_directory() / "hcp_mg.vasp";
    std::ifstream hcp_read_prim_stream(hcp_read_prim_path.string());
    xtal::BasicStructure hcp_read_prim_structure = xtal::BasicStructure::from_poscar_stream(hcp_read_prim_stream);
    hcp_read_prim_ptr.reset(new xtal::BasicStructure(std::move(hcp_read_prim_structure)));

    hcp_3stack_factor_group = xtal::make_factor_group(*hcp_3stack_ptr);
    hcp_primitive_factor_group = xtal::make_factor_group(*hcp_primitive_ptr);
    hcp_read_prim_factor_group = xtal::make_factor_group(*hcp_read_prim_ptr);
    return;
  }

  std::unique_ptr<xtal::BasicStructure> hcp_3stack_ptr;
  std::unique_ptr<xtal::BasicStructure> hcp_primitive_ptr;
  std::unique_ptr<xtal::BasicStructure> hcp_read_prim_ptr;;

  xtal::SymOpVector hcp_3stack_factor_group;
  xtal::SymOpVector hcp_primitive_factor_group;
  xtal::SymOpVector hcp_read_prim_factor_group;;
};

TEST_F(HexagonalSuperStructureTest, MakePrimitive) {

  EXPECT_EQ(hcp_3stack_ptr->basis().size() / 3, hcp_primitive_ptr->basis().size());
}

TEST_F(HexagonalSuperStructureTest, PrimitiveFactorGroup) {
  EXPECT_EQ(hcp_primitive_factor_group.size(), 24);
}

TEST_F(HexagonalSuperStructureTest, ReadPrimitiveFactorGroup) {
  EXPECT_EQ(hcp_read_prim_factor_group.size(), hcp_primitive_factor_group.size());
}

TEST_F(HexagonalSuperStructureTest, SuperStructureFactorGroup) {
  int num_prim_units = hcp_3stack_ptr->basis().size() / hcp_primitive_ptr->basis().size(); //=3
  EXPECT_EQ(hcp_primitive_factor_group.size() * num_prim_units, hcp_3stack_factor_group.size());
}
