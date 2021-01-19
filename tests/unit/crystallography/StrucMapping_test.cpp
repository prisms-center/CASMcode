#include "gtest/gtest.h"

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include "casm/clex/SimpleStructureTools.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/SimpleStrucMapCalculator.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/StrucMapping.hh"

using namespace CASM;

/// \brief Test to ensure that the symmetry invariant mapping costs are
/// identically zero for all values of c/a in an hcp structure
TEST(StrucMapping_Test, Test1) {
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
