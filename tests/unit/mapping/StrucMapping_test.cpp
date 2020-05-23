#include "gtest/gtest.h"
#include <vector>
#include "autotools.hh"

/// What is being tested:
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/StrucMapping.hh"

/// What is being used to test it:
#include "crystallography/TestStructures.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/SimpleStrucMapCalculator.hh"

using namespace CASM;

void print_mapping_nodes(std::set<xtal::MappingNode> const &set) {
  int i = 0;
  for(auto const &el : set) {
    std::cout << "ELEMENT " << ++i << ":\n";
    std::cout << "   cost: " << el.cost << "  bcost: " << el.atomic_node.cost << "  lcost: " << el.lattice_node.cost << "\n"
              << "   translation: " << el.atomic_node.translation.transpose() << "\n"
              << "   isometry: \n" << el.lattice_node.isometry << "\n"
              << "   stretch: \n" << el.lattice_node.stretch << "\n"
              << "   parent: \n" << el.lattice_node.parent.superlattice().lat_column_mat() << "\n"
              << "   cost_mat: \n" << el.atomic_node.cost_mat << "\n"
              << "   partitioned: " << el.is_partitioned << "\n"
              << "   forced_on: \n";
    for(auto const &pr : el.atomic_node.forced_on)
      std::cout << "     (" << pr.first << ", " << pr.second << ")\n";
    std::cout << "   irow: " << el.atomic_node.irow << "\n"
              << "   icol: " << el.atomic_node.icol << "\n"
              << "   assignment: " << el.atomic_node.assignment << "\n"
              << "   displacement: \n" << el.atom_displacement << "\n"
              << "   tot assignment: " << el.atom_permutation << "\n\n-----\n\n";
  }
}
// Generate cubic cell with lat param a and two atoms of species "A" separated by d along [111]
// when d=sqrt(3)a/2, describes BCC
xtal::SimpleStructure map_struc1(double a, double d) {
  xtal::SimpleStructure result;
  result.lat_column_mat.setIdentity();
  result.lat_column_mat *= a;
  result.atom_info.resize(2);
  result.atom_info.names[0] = "A";
  result.atom_info.names[1] = "A";
  result.atom_info.cart_coord(0).setZero();
  result.atom_info.cart_coord(1).setConstant(d / sqrt(3.));
  return result;
}

void k_best_mapping_test(xtal::SimpleStructure const &sstruc, double d) {

  // Store result as factor group of structure
  xtal::SymOpVector fgroup;
  {
    std::string comment("Check for perfect mappings using the best-0 calling convention, without symmetry and with a positive min_cost");
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(sstruc, xtal::Lattice(sstruc.lat_column_mat), 0, xtal::StrucMapping::big_inf(), 1e-3);

    EXPECT_EQ(sym_set.size(), 12) << comment;
    fgroup = adapter::Adapter<xtal::SymOpVector, decltype(sym_set)>()(sym_set);

    //std::cout << "BASE MAPPINGS:\n";
    print_mapping_nodes(sym_set);
  }

  {
    std::string comment("Check for best all mappings better than the pure swap mapping, which has a cost of d^2. There are 8");
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc, fgroup)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(sstruc, xtal::Lattice(sstruc.lat_column_mat), 200, d * d + 1e-6, -1e-3);

    EXPECT_EQ(sym_set.size(), 8) << comment;

    //std::cout << "SUB MAPPINGS:\n";
    //print_mapping_nodes(sym_set);

    EXPECT_NEAR(sym_set.begin()->cost, 0, 1e-6) << comment;

    EXPECT_NEAR(sym_set.rbegin()->cost, 0.5 * d * d, 1e-6) << comment;
  }

}

void sym_mapping_test(xtal::BasicStructure struc, Index N) {
  xtal::SimpleStructure sstruc = xtal::make_simple_structure(struc);

  for(std::string &sp : sstruc.mol_info.names) {
    if(sp == "Va") {
      sp = "A";
    }
  }

  for(std::string &sp : sstruc.atom_info.names) {
    if(sp == "Va") {
      sp = "A";
    }
  }

  //std::cout << "species:  ";
  //for(Index i=0; i<sstruc.atom_info.size(); ++i){
  //std::cout << sstruc.atom_info.coords.col(i).transpose() << "  " << sstruc.atom_info.names[i] << "\n";
  //}
  //std::cout << "\n";

  Eigen::Matrix3i T;
  T.setIdentity();
  T *= 2;

  xtal::SimpleStructure sstruc2 = make_superstructure(T, sstruc);

  {
    std::string comment("Check that we find 8 perfect mapping for a Vol8 non-primitive structure");
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc2)));
    xtal::LatticeNode tnode((xtal::Lattice(sstruc2.lat_column_mat)),
                            (xtal::Lattice(sstruc2.lat_column_mat)),
                            (xtal::Lattice(sstruc2.lat_column_mat)),
                            (xtal::Lattice(sstruc2.lat_column_mat)),
                            sstruc2.atom_info.size());

    auto trans_set = mapper.map_deformed_struc_impose_lattice_node(sstruc2, tnode, 0, xtal::StrucMapping::big_inf(), 1e-3);
    EXPECT_EQ(trans_set.size(), 8) << comment;
  }

  {
    std::string comment("Check for perfect mappings using the best-1 calling convention, without symmetry");
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(sstruc, xtal::Lattice(sstruc.lat_column_mat), 1, xtal::StrucMapping::big_inf(), -1e-3);
    EXPECT_EQ(sym_set.size(), N) << comment;
  }

  {
    std::string comment("Check for perfect mappings using the best-1000 calling convention, without symmetry");
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(sstruc, xtal::Lattice(sstruc.lat_column_mat), 1000, 1e-3, -1e-3);
    EXPECT_EQ(sym_set.size(), N) << comment;
  }

  // Store result as factor group of structure
  xtal::SymOpVector fgroup;
  {
    std::string comment("Check for perfect mappings using the best-0 calling convention, without symmetry and with a positive min_cost");
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(sstruc, xtal::Lattice(sstruc.lat_column_mat), 0, xtal::StrucMapping::big_inf(), 1e-3);

    EXPECT_EQ(sym_set.size(), N) << comment;
    fgroup = adapter::Adapter<xtal::SymOpVector, decltype(sym_set)>()(sym_set);
  }

  {
    std::string comment("Check for perfect mappings of primitive structure onto itself, using symmetry reduction of factor group from previous step.");
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc, fgroup)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(sstruc, xtal::Lattice(sstruc.lat_column_mat), 0, xtal::StrucMapping::big_inf(), 1e-3);

    EXPECT_EQ(sym_set.size(), 1) << comment;
  }

  {
    std::string comment("Check for perfect mappings of non-primitive structure onto primitive, using symmetry reduction of factor group from previous step.");
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc, fgroup)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(sstruc2, xtal::Lattice(sstruc2.lat_column_mat), 0, xtal::StrucMapping::big_inf(), 1e-3);
    EXPECT_EQ(sym_set.size(), 1) << comment;
  }

  {
    std::string comment("Check for perfect mappings of vol-8 non-primitive structure onto itself, using symmetry reduction of factor group from previous step.");
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc2, fgroup)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(sstruc2, xtal::Lattice(sstruc2.lat_column_mat), 0, xtal::StrucMapping::big_inf(), 1e-3);
    EXPECT_EQ(sym_set.size(), 8) << comment;
  }



}

xtal::BasicStructure select_first_allowed_occupant_for_basis(const xtal::BasicStructure &ambiguous_structure) {
  xtal::BasicStructure struc_in_state = ambiguous_structure;
  std::vector<xtal::Site> new_basis;
  for(const xtal::Site &s : struc_in_state.basis()) {
    new_basis.emplace_back(s, s.allowed_occupants()[0]);
  }

  struc_in_state.set_basis(new_basis);
  return struc_in_state;
}


TEST(SymMappingTest1, FCCTernaryPrim) {
  // Read in test PRIM and run tests
  xtal::BasicStructure prim = test::FCC_ternary_prim();
  sym_mapping_test(::select_first_allowed_occupant_for_basis(prim), 48);
}

TEST(SymMappingTest2, ZrOPrim) {
  // Read in test PRIM and run tests
  xtal::BasicStructure prim = test::ZrO_prim();
  sym_mapping_test(::select_first_allowed_occupant_for_basis(prim), 24);
}

TEST(KBestMappingTest, Struc1) {
  // Read in test PRIM and run tests
  k_best_mapping_test(map_struc1(5., 0.5), 0.5);
}


