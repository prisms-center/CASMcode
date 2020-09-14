#include "gtest/gtest.h"
#include "autotools.hh"
#include "Common.hh"

#include "casm/app/ProjectBuilder.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/database/ScelDatabaseTools.hh"
#include "casm/clex/ConfigEnumAllOccupations_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/Superlattice.hh"

#include "crystallography/TestStructures.hh" // for test::ZrO_prim

// ConfigEnumInput
// ---------------
//
// All Configuration enumerators act on a ConfigEnumInput object which specifies a starting
// Configuration used to indicate
// and enumerator specific which provides in CASM are classes that provide iterators which when incremented iteratively
// constuct new objects, typically Supercell or Configuration. When used via the casm command line
// program subcommand `casm enum`, the constructed objects are added to a database for future use.
// When used in C++ code, the constructed objects can be stored in the database or the used in other
// ways.
//
// This example demonstrates enumerating Supercell. There are three related Supercell enumerators:
// - `ScelEnumByProps`: Enumerate Supercell by enumerating superlattices as specifying parameters
//   (CASM::xtal::ScelEnumProps) such as the beginning volume, ending volume, what the unit
//   lattice is (in terms of the prim lattice), and which lattice vectors to enumerate over. This
//   is similar to the example 002_crystallography_superlattice_test.cpp.
// - `ScelEnumByName`: Iterates over Supercell that already exist in the Supercell database by
//   specifying a list of Supercell by name. This is mostly useful as an input to other methods
//   specifying which Supercells to use as input.
// - `ScelEnum`: This enumerator is primarily intended for command line program use and allows use
//   of either `ScelEnumByProps` or `ScelEnumByName` depending on which parameters are passed.
//

// This test fixture class constructs a CASM PrimClex and Supercells for enumeration examples
class ExampleEnumerationZrOConfigEnumAllOccupations : public testing::Test {
protected:

  std::string title;
  std::shared_ptr<CASM::Structure const> shared_prim;
  CASM::ProjectSettings project_settings;
  CASM::PrimClex primclex;

  ExampleEnumerationZrOConfigEnumAllOccupations():
    title("ExampleEnumerationZrOConfigEnumAllOccupations"),
    shared_prim(std::make_shared<CASM::Structure const>(test::ZrO_prim())),
    project_settings(make_default_project_settings(*shared_prim, title)),
    primclex(project_settings, shared_prim) {

    int begin_volume {1};
    int end_volume {5};
    std::string dirs {"abc"};
    Eigen::Matrix3i generating_matrix {Eigen::Matrix3i::Identity()};
    CASM::xtal::ScelEnumProps enumeration_params {begin_volume, end_volume, dirs, generating_matrix};
    bool existing_only = false;

    // The ScelEnumByProps variant that accepts a PrimClex in the constructor inserts Supercells into
    // the Supercell database available at `primclex.db<Supercell>()` as it constructs them.
    CASM::ScelEnumByProps enumerator {primclex, enumeration_params, existing_only};

    // Increments the enumerator iterators to construct all Supercell
    int count = std::distance(enumerator.begin(), enumerator.end());
    EXPECT_EQ(count, 20);
    EXPECT_EQ(primclex.db<Supercell>().size(), 20);
  }

};

TEST_F(ExampleEnumerationZrOConfigEnumAllOccupations, Example1) {
  // Example enumerating all symmetrically unique occupational orderings in each of the supercells
  // in the Supercell database
  std::vector<CASM::Configuration> configurations;
  for(auto const &scel : primclex.db<Supercell>()) {

    // Constructing ConfigEnumInput with a Supercell does:
    CASM::ConfigEnumInput input {scel};
    CASM::ConfigEnumAllOccupations enumerator {input};
    std::copy(enumerator.begin(), enumerator.end(), std::back_inserter(configurations));
  }
  EXPECT_EQ(configurations.size(), 336);
}

TEST_F(ExampleEnumerationZrOConfigEnumAllOccupations, Example2) {

  // Enumerate point and pair cluster perturbations of an ordered structure

  // Configuration with composition Zr_{6}O_{1}

  Eigen::Matrix3l T;
  T << 2, -1, 1,
  1, 1, 1,
  0, 0, 1;

  CASM::Supercell const &supercell = *make_canonical_and_insert(shared_prim, T, primclex.db<Supercell>()).first;

  std::cout << "super lattice: \n" << supercell.lattice().lat_column_mat() << std::endl;
  std::cout << "T:\n" << supercell.sym_info().transformation_matrix_to_super() << std::endl;
}
