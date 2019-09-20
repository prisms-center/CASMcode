#include "gtest/gtest.h"

/// What is being tested:
#include "casm/crystallography/Molecule.hh"

/// What is being used to test it:
#include "casm/misc/CASM_Eigen_math.hh"

using namespace CASM;

TEST(MoleculeTest, AtomPositionTest1) {

  Eigen::Vector3d vec(0.0, 0.2, 0.4);
  double tol(1e-5);

  AtomPosition atom_pos(vec, "A");

  EXPECT_EQ(atom_pos.name(), "A");
  EXPECT_EQ(almost_equal(atom_pos.cart(), vec, tol), true);

}

TEST(MoleculteTest, MoleculeTest1) {

  Eigen::Vector3d vec(0.0, 0.2, 0.4);
  double tol(1e-5);
  bool divisible(true);

  Molecule mol_a = Molecule::make_atom("A");
  Molecule mol_h2o = Molecule("H2O", {AtomPosition(vec, "H"), AtomPosition(vec, "H"), AtomPosition(vec, "O")});
  Molecule mol_Zr2 = Molecule("Zr2", {AtomPosition(vec, "Zr"), AtomPosition(-1.0 * vec, "Zr")}, divisible);
  Molecule mol_va = Molecule::make_vacancy();
  Molecule mol_va2 = Molecule::make_atom("Va"); // will make vacancy with no atoms

  EXPECT_EQ(mol_a.size(), 1);
  EXPECT_EQ(mol_a.name(), "A");
  EXPECT_EQ(mol_a.is_vacancy(), false);
  EXPECT_EQ(mol_a.atom(0).name() == "A", true);

  EXPECT_EQ(mol_h2o.size(), 3);
  EXPECT_EQ(mol_h2o.name(), "H2O");
  EXPECT_EQ(mol_h2o.is_vacancy(), false);
  EXPECT_EQ(mol_h2o.atoms().size(), 3);
  EXPECT_EQ(mol_h2o.atoms().at(0).name(), "H");
  EXPECT_EQ(mol_h2o.atoms().at(1).name(), "H");
  EXPECT_EQ(mol_h2o.atoms().at(2).name(), "O");
  EXPECT_EQ(mol_h2o.atom(0).name(), "H");
  EXPECT_EQ(mol_h2o.atom(1).name(), "H");
  EXPECT_EQ(mol_h2o.atom(2).name(), "O");
  EXPECT_EQ(mol_h2o.is_divisible(), false);
  EXPECT_EQ(mol_h2o.is_indivisible(), true);

  EXPECT_EQ(mol_Zr2.size(), 2);
  EXPECT_EQ(mol_Zr2.name(), "Zr2");
  EXPECT_EQ(mol_Zr2.is_vacancy(), false);
  EXPECT_EQ(mol_Zr2.is_divisible(), true);
  EXPECT_EQ(mol_Zr2.is_indivisible(), false);

  //EXPECT_EQ(mol_va.size(), 0); // if Molecule empty for Va
  EXPECT_EQ(mol_va.size(), 1);
  EXPECT_EQ(mol_va.name(), "Va");
  EXPECT_EQ(mol_va.is_vacancy(), true);

  //EXPECT_EQ(mol_va2.size(), 0); // if Molecule empty for Va
  EXPECT_EQ(mol_va2.size(), 1);
  EXPECT_EQ(mol_va2.name(), "Va");
  EXPECT_EQ(mol_va2.is_vacancy(), true);
}

