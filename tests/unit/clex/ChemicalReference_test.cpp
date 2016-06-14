#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clex/ChemicalReference.hh"

/// What is being used to test it:

#include "casm/app/ProjectBuilder.hh"
#include "Common.hh"
#include "FCCTernaryProj.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(ChemicalReferenceTest)

BOOST_AUTO_TEST_CASE(Basics) {

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
  ChemicalReference chem_ref(
    prim,
    ref_state.begin(),
    ref_state.end(),
    tol
  );

  // check the global reference hyperplane
  //std::cout << "global ref: " << chem_ref.global().transpose() << std::endl;
  BOOST_CHECK_EQUAL(almost_equal(chem_ref.global(), Eigen::Vector3d(-1.0, -2.0, -4.0), tol), true);

  // this should not be allowed to compile:
  //chem_ref.global() = Eigen::Vector3d::Zero();

  // try writing to json
  jsonParser json;
  to_json(chem_ref, json);
  BOOST_CHECK_EQUAL(json.contains("global"), true);

  // try reading from json
  ChemicalReference chem_ref2 = json.get<ChemicalReference>(prim, tol);
  BOOST_CHECK_EQUAL(almost_equal(chem_ref2.global(), Eigen::Vector3d(-1.0, -2.0, -4.0), tol), true);


  // try writing to json, as in project file
  write_chemical_reference(chem_ref, json);
  BOOST_CHECK_EQUAL(json.contains("chemical_reference"), true);

  // try reading from json, as in a project file
  ChemicalReference chem_ref3 = read_chemical_reference(json, prim, tol);
  BOOST_CHECK_EQUAL(almost_equal(chem_ref3.global(), Eigen::Vector3d(-1.0, -2.0, -4.0), tol), true);


  // check the print formatting
  std::stringstream ss;
  ChemicalReferencePrinter p(ss, chem_ref);
  p.print_all();
  //std::cout << ss.str() << std::endl;

}

BOOST_AUTO_TEST_CASE(HyperplaneConstructor) {

  Structure prim(test::FCC_ternary_prim());

  // construct the ChemicalReference, providing only a global reference
  ChemicalReference chem_ref(
    prim,
    Eigen::Vector3d(-1.0, -2.0, -4.0)
  );

  // check the global reference hyperplane
  //std::cout << "global ref: " << chem_ref.global().transpose() << std::endl;
  BOOST_CHECK_EQUAL(almost_equal(chem_ref.global(), Eigen::Vector3d(-1.0, -2.0, -4.0), TOL), true);

  // this should not be allowed to compile:
  //chem_ref.global() = Eigen::Vector3d::Zero();

  // try writing to json
  jsonParser json;
  to_json(chem_ref, json);
  BOOST_CHECK_EQUAL(json.contains("global"), true);
  //std::cout << json << std::endl;

  // try writing to json, as in project file
  write_chemical_reference(chem_ref, json);
  BOOST_CHECK_EQUAL(json.contains("chemical_reference"), true);
  //std::cout << json << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
