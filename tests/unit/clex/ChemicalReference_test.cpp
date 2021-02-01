#include "gtest/gtest.h"

/// What is being tested:
#include "casm/clex/ChemicalReference_impl.hh"
#include "casm/clex/io/json/ChemicalReference_json_io.hh"

/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/app/ProjectBuilder.hh"
#include "casm/crystallography/Structure.hh"

using namespace CASM;

TEST(ChemicalReferenceTest, Basics) {
  Structure prim(test::FCC_ternary_prim());
  double tol = TOL;

  // initialize reference states, just using energies for the end members
  std::vector<ChemicalReferenceState> ref_state(3);
  ref_state[0].species_num["A"] = 1.0;
  ref_state[0].energy_per_species = -1.0;

  ref_state[1].species_num["B"] = 1.0;
  ref_state[1].energy_per_species = -2.0;

  ref_state[2].species_num["C"] = 1.0;
  ref_state[2].energy_per_species = -4.0;

  // construct the ChemicalReference, providing only a global reference
  ChemicalReference chem_ref(prim, ref_state.begin(), ref_state.end(), tol);

  // check the global reference hyperplane
  // std::cout << "global ref: " << chem_ref.global().transpose() << std::endl;
  EXPECT_EQ(
      almost_equal(chem_ref.global(), Eigen::Vector3d(-1.0, -2.0, -4.0), tol),
      true);

  // this should not be allowed to compile:
  // chem_ref.global() = Eigen::Vector3d::Zero();

  // try writing to json
  jsonParser json;
  to_json(chem_ref, json);
  EXPECT_EQ(json.contains("global"), true);

  // try reading from json
  ChemicalReference chem_ref2 = json.get<ChemicalReference>(prim, tol);
  EXPECT_EQ(
      almost_equal(chem_ref2.global(), Eigen::Vector3d(-1.0, -2.0, -4.0), tol),
      true);

  // try writing to json, as in project file
  to_json(chem_ref, json["chemical_reference"]);
  EXPECT_EQ(json.contains("chemical_reference"), true);

  // try reading from json, as in a project file
  ChemicalReference chem_ref3 =
      json["chemical_reference"].get<ChemicalReference>(prim, tol);

  EXPECT_EQ(
      almost_equal(chem_ref3.global(), Eigen::Vector3d(-1.0, -2.0, -4.0), tol),
      true);

  // check the print formatting
  std::stringstream ss;
  ChemicalReferencePrinter p(ss, chem_ref);
  p.print_all();
  // std::cout << ss.str() << std::endl;
}

TEST(ChemicalReferenceTest, HyperplaneConstructor) {
  Structure prim(test::FCC_ternary_prim());

  // construct the ChemicalReference, providing only a global reference
  ChemicalReference chem_ref(prim, Eigen::Vector3d(-1.0, -2.0, -4.0));

  // check the global reference hyperplane
  // std::cout << "global ref: " << chem_ref.global().transpose() << std::endl;
  EXPECT_EQ(
      almost_equal(chem_ref.global(), Eigen::Vector3d(-1.0, -2.0, -4.0), TOL),
      true);

  // this should not be allowed to compile:
  // chem_ref.global() = Eigen::Vector3d::Zero();

  // try writing to json
  jsonParser json;
  to_json(chem_ref, json);
  EXPECT_EQ(json.contains("global"), true);
  // std::cout << json << std::endl;

  // try writing to json, as in project file
  to_json(chem_ref, json["chemical_reference"]);
  EXPECT_EQ(json.contains("chemical_reference"), true);
  // std::cout << json << std::endl;
}
