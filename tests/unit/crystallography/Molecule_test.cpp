#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/crystallography/Molecule.hh"

/// What is being used to test it:
#include "casm/misc/CASM_Eigen_math.hh"

using namespace CASM;


BOOST_AUTO_TEST_SUITE(MoleculeTest)

BOOST_AUTO_TEST_CASE(SpecieTest1) {
  AtomSpecies species_a("A");
  AtomSpecies species_a2("A");
  AtomSpecies species_b("B");

  BOOST_CHECK_EQUAL(species_a.name(), "A");
  BOOST_CHECK_EQUAL(species_a2.name(), "A");
  BOOST_CHECK_EQUAL(species_b.name(), "B");
  BOOST_CHECK_EQUAL(species_a < species_a2, false);
  BOOST_CHECK_EQUAL(species_a == species_a2, true);
  BOOST_CHECK_EQUAL(species_a < species_b, true);
  BOOST_CHECK_EQUAL(species_a > species_b, false);
  BOOST_CHECK_EQUAL(species_a == species_b, false);

}

BOOST_AUTO_TEST_CASE(AtomPositionTest1) {

  Eigen::Vector3d vec(0.0, 0.2, 0.4);
  AtomSpecies species_a("A");
  double tol(1e-5);

  AtomPosition atom_pos_a(vec, species_a);
  AtomPosition atom_pos_a2(vec, "A");
  AtomPosition atom_pos_b(vec, "B");

  BOOST_CHECK_EQUAL(atom_pos_a.name(), "A");
  BOOST_CHECK_EQUAL(atom_pos_a.species() == species_a, true);
  BOOST_CHECK_EQUAL(atom_pos_a2.species() == species_a, true);
  BOOST_CHECK_EQUAL(atom_pos_b.species() != species_a, true);
  BOOST_CHECK_EQUAL(almost_equal(atom_pos_a.cart(), vec, tol), true);

}

BOOST_AUTO_TEST_CASE(MoleculeTest1) {

  Eigen::Vector3d vec(0.0, 0.2, 0.4);
  double tol(1e-5);
  bool divisible(true);

  Molecule mol_a = Molecule::make_atom("A");
  Molecule mol_h2o = Molecule("H2O", {AtomPosition(vec, "H"), AtomPosition(vec, "H"), AtomPosition(vec, "O")});
  Molecule mol_Zr2 = Molecule("Zr2", {AtomPosition(vec, "Zr"), AtomPosition(-1.0 * vec, "Zr")}, divisible);
  Molecule mol_va = Molecule::make_vacancy();
  Molecule mol_va2 = Molecule::make_atom("Va"); // will make vacancy with no atoms

  BOOST_CHECK_EQUAL(mol_a.size(), 1);
  BOOST_CHECK_EQUAL(mol_a.name(), "A");
  BOOST_CHECK_EQUAL(mol_a.is_vacancy(), false);
  BOOST_CHECK_EQUAL(mol_a.atom(0).species() == AtomSpecies("A"), true);

  BOOST_CHECK_EQUAL(mol_h2o.size(), 3);
  BOOST_CHECK_EQUAL(mol_h2o.name(), "H2O");
  BOOST_CHECK_EQUAL(mol_h2o.is_vacancy(), false);
  BOOST_CHECK_EQUAL(mol_h2o.atoms().size(), 3);
  BOOST_CHECK_EQUAL(mol_h2o.atoms().at(0).species() == AtomSpecies("H"), true);
  BOOST_CHECK_EQUAL(mol_h2o.atoms().at(1).species() == AtomSpecies("H"), true);
  BOOST_CHECK_EQUAL(mol_h2o.atoms().at(2).species() == AtomSpecies("O"), true);
  BOOST_CHECK_EQUAL(mol_h2o.atom(0).species() == AtomSpecies("H"), true);
  BOOST_CHECK_EQUAL(mol_h2o.atom(1).species() == AtomSpecies("H"), true);
  BOOST_CHECK_EQUAL(mol_h2o.atom(2).species() == AtomSpecies("O"), true);
  BOOST_CHECK_EQUAL(mol_h2o.is_divisible(), false);
  BOOST_CHECK_EQUAL(mol_h2o.is_indivisible(), true);

  BOOST_CHECK_EQUAL(mol_Zr2.size(), 2);
  BOOST_CHECK_EQUAL(mol_Zr2.name(), "Zr2");
  BOOST_CHECK_EQUAL(mol_Zr2.is_vacancy(), false);
  BOOST_CHECK_EQUAL(mol_Zr2.is_divisible(), true);
  BOOST_CHECK_EQUAL(mol_Zr2.is_indivisible(), false);

  //BOOST_CHECK_EQUAL(mol_va.size(), 0); // if Molecule empty for Va
  BOOST_CHECK_EQUAL(mol_va.size(), 1);
  BOOST_CHECK_EQUAL(mol_va.name(), "Va");
  BOOST_CHECK_EQUAL(mol_va.is_vacancy(), true);

  //BOOST_CHECK_EQUAL(mol_va2.size(), 0); // if Molecule empty for Va
  BOOST_CHECK_EQUAL(mol_va2.size(), 1);
  BOOST_CHECK_EQUAL(mol_va2.name(), "Va");
  BOOST_CHECK_EQUAL(mol_va2.is_vacancy(), true);


}


BOOST_AUTO_TEST_SUITE_END()
