#include "casm/crystallography/StrucMapping.hh"

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include "autotools.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/SimpleStrucMapCalculator.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/io/VaspIO.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "crystallography/TestStructures.hh"
#include "gtest/gtest.h"

using namespace CASM;

void print_mapping_nodes(std::set<xtal::MappingNode> const &set) {
  int i = 0;
  for (auto const &el : set) {
    std::cout << "ELEMENT " << ++i << ":\n";
    std::cout << "   cost: " << el.cost << "  bcost: " << el.atomic_node.cost
              << "  lcost: " << el.lattice_node.cost << "\n"
              << "   translation: " << el.atomic_node.translation.transpose()
              << "\n"
              << "   isometry: \n"
              << el.lattice_node.isometry << "\n"
              << "   stretch: \n"
              << el.lattice_node.stretch << "\n"
              << "   parent: \n"
              << el.lattice_node.parent.superlattice().lat_column_mat() << "\n"
              << "   cost_mat: \n"
              << el.atomic_node.cost_mat << "\n"
              << "   partitioned: " << el.is_partitioned << "\n"
              << "   forced_on: \n";
    for (auto const &pr : el.atomic_node.forced_on)
      std::cout << "     (" << pr.first << ", " << pr.second << ")\n";
    std::cout << "   irow: " << el.atomic_node.irow << "\n"
              << "   icol: " << el.atomic_node.icol << "\n"
              << "   assignment: " << el.atomic_node.assignment << "\n"
              << "   displacement: \n"
              << el.atom_displacement << "\n"
              << "   tot assignment: " << el.atom_permutation
              << "\n\n-----\n\n";
  }
}
// Generate cubic cell with lat param a and two atoms of species "A" separated
// by d along [111] when d=sqrt(3)a/2, describes BCC
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
    std::string comment(
        "Check for perfect mappings using the best-0 calling convention, "
        "without symmetry and with a positive min_cost");
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(
        sstruc, xtal::Lattice(sstruc.lat_column_mat), 0,
        xtal::StrucMapping::big_inf(), 1e-3);

    EXPECT_EQ(sym_set.size(), 12) << comment;
    fgroup = adapter::Adapter<xtal::SymOpVector, decltype(sym_set)>()(sym_set);

    // std::cout << "BASE MAPPINGS:\n";
  }

  {
    std::string comment(
        "Check for all mappings better than the pure swap "
        "mapping, which has a cost of 0.5*d^2. Without "
        "considering symmetry of child structure there are 8.");
    xtal::StrucMapper mapper(xtal::SimpleStrucMapCalculator(sstruc, fgroup),
                             0.5, 0.5, xtal::StrucMapper::robust);
    auto sym_set = mapper.map_deformed_struc_impose_lattice(
        sstruc, xtal::Lattice(sstruc.lat_column_mat), 0,
        xtal::StrucMapping::big_inf(), 0.5 * d * d + 1e-6);

    EXPECT_EQ(sym_set.size(), 8) << comment;

    // std::cout << "SUB MAPPINGS:\n";
    // print_mapping_nodes(sym_set);

    EXPECT_NEAR(sym_set.begin()->cost, 0, 1e-6) << comment;

    EXPECT_NEAR(sym_set.rbegin()->cost, 0.5 * d * d, 1e-6) << comment;
  }

  {
    std::string comment(
        "Check for all mappings better than the pure swap "
        "mapping, which has a cost of 0.5 * d^2. Considering "
        "symmetry of child structure, there are 4.");
    xtal::StrucMapper mapper(xtal::SimpleStrucMapCalculator(sstruc, fgroup),
                             0.5, 0.5, xtal::StrucMapper::robust);
    auto sym_set = mapper.map_deformed_struc_impose_lattice(
        sstruc, xtal::Lattice(sstruc.lat_column_mat), 0,
        xtal::StrucMapping::big_inf(), 0.5 * d * d + 1e-6, false, fgroup);

    EXPECT_EQ(sym_set.size(), 4) << comment;

    // std::cout << "SUB MAPPINGS:\n";
    // print_mapping_nodes(sym_set);

    EXPECT_NEAR(sym_set.begin()->cost, 0, 1e-6) << comment;

    EXPECT_NEAR(sym_set.rbegin()->cost, 0.5 * d * d, 1e-6) << comment;
  }
}

void sym_mapping_test(xtal::BasicStructure struc, Index N) {
  xtal::SimpleStructure sstruc = xtal::make_simple_structure(struc);

  for (std::string &sp : sstruc.mol_info.names) {
    if (sp == "Va") {
      sp = "A";
    }
  }

  for (std::string &sp : sstruc.atom_info.names) {
    if (sp == "Va") {
      sp = "A";
    }
  }

  // std::cout << "species:  ";
  // for(Index i=0; i<sstruc.atom_info.size(); ++i){
  // std::cout << sstruc.atom_info.coords.col(i).transpose() << "  " <<
  // sstruc.atom_info.names[i] << "\n";
  //}
  // std::cout << "\n";

  Eigen::Matrix3i T;
  T.setIdentity();
  T *= 2;

  xtal::SimpleStructure sstruc2 = make_superstructure(T, sstruc);

  {
    std::string comment(
        "Check that we find 8 perfect mapping for a Vol8 "
        "non-primitive structure");
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc2)));
    xtal::LatticeNode tnode((xtal::Lattice(sstruc2.lat_column_mat)),
                            (xtal::Lattice(sstruc2.lat_column_mat)),
                            (xtal::Lattice(sstruc2.lat_column_mat)),
                            (xtal::Lattice(sstruc2.lat_column_mat)),
                            sstruc2.atom_info.size());

    auto trans_set = mapper.map_deformed_struc_impose_lattice_node(
        sstruc2, tnode, 0, xtal::StrucMapping::big_inf(), 1e-3);
    EXPECT_EQ(trans_set.size(), 8) << comment;
  }

  {
    std::string comment(
        "Check for perfect mappings using the best-1 calling "
        "convention, without symmetry");
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(
        sstruc, xtal::Lattice(sstruc.lat_column_mat), 1,
        xtal::StrucMapping::big_inf(), -1e-3);
    EXPECT_EQ(sym_set.size(), N) << comment;
  }

  {
    std::string comment(
        "Check for perfect mappings using the best-1000 "
        "calling convention, without symmetry");
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(
        sstruc, xtal::Lattice(sstruc.lat_column_mat), 1000, 1e-3, -1e-3);
    EXPECT_EQ(sym_set.size(), N) << comment;
  }

  // Store result as factor group of structure
  xtal::SymOpVector fgroup;
  {
    std::string comment(
        "Check for perfect mappings using the best-0 calling convention, "
        "without symmetry and with a positive min_cost");
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(
        sstruc, xtal::Lattice(sstruc.lat_column_mat), 0,
        xtal::StrucMapping::big_inf(), 1e-3);

    EXPECT_EQ(sym_set.size(), N) << comment;
    fgroup = adapter::Adapter<xtal::SymOpVector, decltype(sym_set)>()(sym_set);
  }

  {
    std::string comment(
        "Check for perfect mappings of primitive structure onto itself, using "
        "symmetry reduction of factor group from previous step.");
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc, fgroup)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(
        sstruc, xtal::Lattice(sstruc.lat_column_mat), 0,
        xtal::StrucMapping::big_inf(), 1e-3);

    EXPECT_EQ(sym_set.size(), 1) << comment;
  }

  {
    std::string comment(
        "Check for perfect mappings of non-primitive structure onto primitive, "
        "using symmetry reduction of factor group from previous step.");
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc, fgroup)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(
        sstruc2, xtal::Lattice(sstruc2.lat_column_mat), 0,
        xtal::StrucMapping::big_inf(), 1e-3);
    EXPECT_EQ(sym_set.size(), 1) << comment;
  }

  {
    std::string comment(
        "Check for perfect mappings of vol-8 non-primitive structure onto "
        "itself, using symmetry reduction of factor group from previous step.");
    xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc2, fgroup)));
    auto sym_set = mapper.map_deformed_struc_impose_lattice(
        sstruc2, xtal::Lattice(sstruc2.lat_column_mat), 0,
        xtal::StrucMapping::big_inf(), 1e-3);
    EXPECT_EQ(sym_set.size(), 8) << comment;
  }
}

xtal::BasicStructure select_first_allowed_occupant_for_basis(
    const xtal::BasicStructure &ambiguous_structure) {
  xtal::BasicStructure struc_in_state = ambiguous_structure;
  std::vector<xtal::Site> new_basis;
  for (const xtal::Site &s : struc_in_state.basis()) {
    new_basis.emplace_back(s, s.allowed_occupants()[0]);
  }

  struc_in_state.set_basis(new_basis);
  return struc_in_state;
}

TEST(SymMappingTest, FCCTernaryPrim) {
  // Read in test PRIM and run tests
  xtal::BasicStructure prim = test::FCC_ternary_prim();
  sym_mapping_test(::select_first_allowed_occupant_for_basis(prim), 48);
}

TEST(SymMappingTest, ZrOPrim) {
  // Read in test PRIM and run tests
  xtal::BasicStructure prim = test::ZrO_prim();
  sym_mapping_test(::select_first_allowed_occupant_for_basis(prim), 24);
}

TEST(KBestMappingTest, Struc1) {
  // Read in test PRIM and run tests

  k_best_mapping_test(map_struc1(pow(8. * M_PI / 3, 1. / 3.), 0.3), 0.3);
}

TEST(SelfMappingTest, MultipleAllowedOccupants) {
  Eigen::Matrix3d fcc_lat_vecs;
  fcc_lat_vecs << 0, 2, 2, 0, 2, 0, 2, 2, 0;

  auto make_fcc = [=](const std::string &occupant) {
    xtal::SimpleStructure fcc;
    fcc.lat_column_mat = fcc_lat_vecs;
    fcc.atom_info.resize(1);
    fcc.atom_info.names[0] = occupant;
    return fcc;
  };

  auto parent = make_fcc("Z");
  auto fcc_X = make_fcc("X");
  auto fcc_Y = make_fcc("Y");

  xtal::SimpleStrucMapCalculator iface(
      parent, {CASM::xtal::SymOp::identity()},
      CASM::xtal::SimpleStructure::SpeciesMode::ATOM, {{"Y", "X"}});
  xtal::StrucMapper mapper(iface);

  auto X_results = mapper.map_deformed_struc(fcc_X);
  auto Y_results = mapper.map_deformed_struc(fcc_Y);
  auto Z_results = mapper.map_deformed_struc(parent);

  EXPECT_EQ(Z_results.size(),
            Index{0});  // Z is not an allowed occupant in the mapper
  EXPECT_EQ(X_results.size(), Index{48});
  EXPECT_EQ(X_results.size(), Y_results.size());
}

/// \brief Test to ensure that the symmetry invariant mapping costs are
/// identically zero for all values of c/a in an hcp structure
TEST(SymInvariantMappingTest, Hcp) {
  // Initial testing values and ranges:
  int num_points = 10;
  double ca_min(0.5), ca_max(2.0);
  // std::vector<std::vector<double>> map_cost;
  std::vector<double> map_cost;

  // Initialize the hcp structure
  std::string hcp_struc_string("");
  hcp_struc_string += "Ti2\n";
  hcp_struc_string += "1.0\n";
  hcp_struc_string += "2.871000 -0.000000 0.000000\n";
  hcp_struc_string += "-1.435500 2.486359 0.000000\n";
  hcp_struc_string += "0.000000 0.000000 4.635000\n";
  hcp_struc_string += "Ti\n";
  hcp_struc_string += "2\n";
  hcp_struc_string += "direct\n";
  hcp_struc_string += "0.666667 0.333333 0.250000 Ti\n";
  hcp_struc_string += "0.333333 0.666667 0.750000 Ti\n";

  // Convert the structure to a SimpleStructure and a BasicStructure
  std::stringstream hcp_string_stream;
  hcp_string_stream << hcp_struc_string;
  xtal::BasicStructure basic_struc =
      xtal::BasicStructure::from_poscar_stream(hcp_string_stream);
  xtal::SimpleStructure simple_struc = xtal::make_simple_structure(basic_struc);

  // Get the factor group of the structure
  auto parent_fg = xtal::make_factor_group(basic_struc, 0.0001);
  EXPECT_EQ(parent_fg.size(), 24)
      << "The obtained factor group for a hcp crystal structure had "
      << parent_fg.size() << " operations. 24 operations are expected";

  // Initialize the defaults for structure mapping
  double lattice_weight(0.5), max_vol_change(2), cost_tol(0.00001),
      min_va_frac(0.0), max_va_frac(0.0);
  int options(xtal::StrucMapper::robust);

  // Initialize a symmetrized  strucmapper
  xtal::StrucMapper sym_struc_map(
      xtal::SimpleStrucMapCalculator(simple_struc, parent_fg,
                                     xtal::SimpleStructure::SpeciesMode::ATOM,
                                     allowed_molecule_names(basic_struc)),
      lattice_weight, max_vol_change, options, cost_tol, min_va_frac,
      max_va_frac);
  sym_struc_map.set_symmetrize_lattice_cost(true);

  double ca_inc = (ca_max - ca_min) / num_points;
  Index k_best = 1;
  double max_cost(1e10), min_cost(-1);
  Eigen::Vector3d parent_lattice_parameters =
      simple_struc.lat_column_mat.colwise().norm();
  bool use_child_sym = true;

  // Loop through all c/a ratios and calculate the mapping cost
  for (double ca = ca_min; ca <= ca_max;) {
    // Strain the structure and prepare it for the mapping routines
    xtal::SimpleStructure strained_struc =
        xtal::make_simple_structure(basic_struc);
    Eigen::Matrix3d _tmp_strained_lattice, _tmp_deformation_tensor;
    _tmp_strained_lattice = strained_struc.lat_column_mat;
    _tmp_strained_lattice.col(2) =
        (ca * parent_lattice_parameters[0] / parent_lattice_parameters[2]) *
        simple_struc.lat_column_mat.col(2);
    _tmp_deformation_tensor =
        _tmp_strained_lattice * simple_struc.lat_column_mat.inverse();
    strained_struc.deform_coords(_tmp_deformation_tensor);
    auto child_fg = xtal::make_factor_group(basic_struc);

    // Calculate the symmetric map cost
    auto sym_tresult = sym_struc_map.map_deformed_struc(
        strained_struc, k_best, max_cost, min_cost, false,
        use_child_sym ? child_fg : decltype(child_fg){child_fg[0]});

    // There should only be a single viable map
    EXPECT_EQ(sym_tresult.size(), 1)
        << "Expected only a single viable map, however " << sym_tresult.size()
        << " maps were found";

    for (auto result_itr = sym_tresult.begin(); result_itr != sym_tresult.end();
         ++result_itr) {
      map_cost.push_back(result_itr->cost);
    }
    ca += ca_inc;
  }

  // Ensure that the range of maps are within a numerical tolerance
  double max_map_cost = *(std::max_element(map_cost.begin(), map_cost.end()));
  double min_map_cost = *(std::min_element(map_cost.begin(), map_cost.end()));
  double MAP_TOL = 0.0001;
  EXPECT_LE(std::abs(max_map_cost), MAP_TOL)
      << "Expected the max mapping cost to be less than " << MAP_TOL
      << ". The maximum mapping cost obtained was " << max_map_cost;
  EXPECT_LE(std::abs(min_map_cost), MAP_TOL)
      << "Expected the min mapping cost to be less than " << MAP_TOL
      << ". The minimum mapping cost obtained was " << max_map_cost;
}

/// \brief Test to ensure that the symmetry invariant mapping costs are
/// identically zero for all values of c/a in an hcp structure
TEST(SymInvariantMappingTest, shuffle) {
  // Setup the initial structure
  std::string orthorhombic_string("");
  orthorhombic_string += "title\n";
  orthorhombic_string += "1.0\n";
  orthorhombic_string += "0.000000000000 1.643000006700 -2.323552846900\n";
  orthorhombic_string += "0.000000000000 1.643000006700  2.323552846900\n";
  orthorhombic_string += "4.647105693800 0.000000000000  0.000000000000\n";
  orthorhombic_string += "H\n";
  orthorhombic_string += "2\n";
  orthorhombic_string += "direct\n";
  orthorhombic_string += "0.000000000000 0.000000000000 0.000000000000\n";
  orthorhombic_string += "0.416666985000 0.583333015000 0.500000000000\n";

  // Convert the structure to a SimpleStructure and a BasicStructure
  std::stringstream ortho_string_stream;
  ortho_string_stream << orthorhombic_string;
  xtal::BasicStructure basic_struc =
      xtal::BasicStructure::from_poscar_stream(ortho_string_stream);
  xtal::SimpleStructure simple_struc = xtal::make_simple_structure(basic_struc);

  // Initial testing values and ranges:
  int num_points = 20;
  double max_shuffle_amplitude(-1.0);
  Eigen::Vector3d unit_shuffle =
      simple_struc.lat_column_mat.col(1) - simple_struc.lat_column_mat.col(0);
  unit_shuffle = unit_shuffle / unit_shuffle.norm();
  Index shuffle_atom_idx = 1;
  double shuffle_inc = max_shuffle_amplitude / num_points;
  std::vector<double> map_cost;

  // Get the factor group of the structure
  auto parent_fg = xtal::make_factor_group(basic_struc, 0.0001);
  EXPECT_EQ(parent_fg.size(), 8)
      << "Expected the factor group size to be 8, instead " << parent_fg.size()
      << " operations were found" << std::endl;

  // Initialize the defaults for structure mapping
  double lattice_weight(0.5), max_vol_change(2), cost_tol(0.00001),
      min_va_frac(0.0), max_va_frac(0.0);
  Index k_best = 1;
  double max_cost(1e10), min_cost(-1);
  bool use_child_sym = true;
  int options(xtal::StrucMapper::robust);

  // Initialize a symmetrized  strucmapper
  xtal::StrucMapper sym_struc_map(
      xtal::SimpleStrucMapCalculator(simple_struc, parent_fg,
                                     xtal::SimpleStructure::SpeciesMode::ATOM,
                                     allowed_molecule_names(basic_struc)),
      lattice_weight, max_vol_change, options, cost_tol, min_va_frac,
      max_va_frac);
  sym_struc_map.set_symmetrize_atomic_cost(true);

  auto child_fg = xtal::make_factor_group(basic_struc);

  // Loop through all shuffles and calculate the mapping cost
  for (double shuffle = 0.0; shuffle >= max_shuffle_amplitude;) {
    // Strain the structure and prepare it for the mapping routines
    xtal::SimpleStructure shuffled_struc =
        xtal::make_simple_structure(basic_struc);
    shuffled_struc.atom_info.cart_coord(shuffle_atom_idx) =
        shuffled_struc.atom_info.cart_coord(shuffle_atom_idx) +
        shuffle * unit_shuffle;

    // Calculate the symmetric map cost
    auto sym_tresult = sym_struc_map.map_deformed_struc(
        shuffled_struc, k_best, max_cost, min_cost, false,
        use_child_sym ? child_fg : decltype(child_fg){child_fg[0]});

    for (auto result_itr = sym_tresult.begin(); result_itr != sym_tresult.end();
         ++result_itr) {
      map_cost.push_back(result_itr->atomic_node.cost);
    }
    shuffle += shuffle_inc;
  }

  // Ensure that the range of maps are within a numerical tolerance
  double max_map_cost = *(std::max_element(map_cost.begin(), map_cost.end()));
  double min_map_cost = *(std::min_element(map_cost.begin(), map_cost.end()));
  double MAP_TOL = 0.00001;
  EXPECT_LE(std::abs(max_map_cost), MAP_TOL)
      << "Expected the max mapping cost to be less than " << MAP_TOL
      << ". The maximum mapping cost obtained was " << max_map_cost;
  EXPECT_LE(std::abs(min_map_cost), MAP_TOL)
      << "Expected the min mapping cost to be less than " << MAP_TOL
      << ". The minimum mapping cost obtained was " << max_map_cost;
}
