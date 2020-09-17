#include "gtest/gtest.h"
#include "autotools.hh"
#include "casm/app/AppIO_impl.hh" // for orbit printing
#include "casm/casm_io/container/json_io.hh"
#include "casm/clex/FillSupercell.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "casm/clusterography/SupercellClusterOrbits_impl.hh"
#include "casm/clusterography/io/stream/IntegralCluster_stream_io.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/enumerator/ClusterSitesSelector_impl.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/symmetry/SubOrbits_impl.hh"
#include "crystallography/TestStructures.hh"

using namespace CASM;
using namespace test;

// helper functions
namespace {

  /// Make transformation_matrix_to_super for conventional 4 atom FCC cell, from the 1 atom primitive
  Eigen::Matrix3i _fcc_conventional_transf_mat();

  /// Make L12 (A3B1), in conventional 4 atom FCC cell
  Configuration _make_configuration_L12(std::shared_ptr<CASM::Structure const> shared_prim);

  /// 3x3x3 super configuration of configuration_L12
  Configuration _make_configuration_L12_3x3x3(Configuration const &_configuration_L12);

  /// NxNxN super configuration of configuration_L12
  Configuration _make_configuration_L12_big(long N, Configuration const &_configuration_L12);

  /// Check how many ConfigEnumInput have a matching number of each type of occupant on the selected sites
  Index count_selected(std::vector<ConfigEnumInput> const &with_cluster_sites, std::vector<Index> expected_occ_count);

}

class ClusterSitesSelectorTest : public testing::Test {
protected:

  std::shared_ptr<CASM::Structure const> shared_prim;
  Configuration configuration_L12;        // L12 (A3B1), in conventional 4 atom FCC cell
  Configuration configuration_L12_3x3x3;  // 3x3x3 superconfiguration of configuration_L12
  Configuration configuration_L12_big;    // big superconfiguration of configuration_L12

  ClusterSitesSelectorTest():
    shared_prim(std::make_shared<CASM::Structure const>(test::FCC_ternary_prim())),
    configuration_L12(_make_configuration_L12(shared_prim)),
    configuration_L12_3x3x3(_make_configuration_L12_3x3x3(configuration_L12)),
    configuration_L12_big(_make_configuration_L12_big(3, configuration_L12)) {
    EXPECT_EQ(configuration_L12.size(), 4);
    EXPECT_EQ(configuration_L12_3x3x3.size(), 108);
  }

};

TEST_F(ClusterSitesSelectorTest, Test1) {

  // Example, generating Configuration with unique cluster sites selected:
  // (starting from "periodic_max_length" clusters)
  // - L12 (A3B1) in conventional FCC supercell (4 atoms)

  ConfigEnumInput reference_config_enum_input {configuration_L12, {}};

  // double a = shared_prim->lattice()[0].norm();
  double a = configuration_L12.supercell().lattice()[0].norm();

  std::vector<double> max_length {
    0, // null cluster
    0, // point clusters
    3 * a + TOL, // a * sqrt(2.0) / 2.0 + TOL // 1NN pair clusters
    3 * a + TOL // a * sqrt(2.0) / 2.0 + TOL // 1NN pair clusters
  };

  // TODO: update this

  PeriodicMaxLengthClusterSpecs cluster_specs {
    shared_prim,
    shared_prim->factor_group(),
    alloy_sites_filter,
    max_length}; // max length, for point clusters only (null orbit branch, point orbit branch)
  auto orbits = cluster_specs.make_periodic_orbits(log());

  SupercellSymInfo const &supercell_sym_info = configuration_L12.supercell().sym_info();
  std::vector<IntegralCluster> scel_generators;
  make_suborbit_generators(
    shared_prim->factor_group().begin(),
    shared_prim->factor_group().end(),
    supercell_sym_info.factor_group().begin(),
    supercell_sym_info.factor_group().end(),
    prototype_iterator(orbits.begin()),
    prototype_iterator(orbits.end()),
    cluster_specs.sym_compare,
    std::back_inserter(scel_generators));

  // Construct OrbitPrinterOptions:
  OrbitPrinterOptions printer_options;
  printer_options.coord_type = CART;
  printer_options.delim = char();
  printer_options.orbit_print_mode = ORBIT_PRINT_MODE::PROTO;
  Printer<IntegralCluster> printer {printer_options};

  Index i;
  i = 0;
  log() << "supercell_sym_info.factor_group().size(): "
        << supercell_sym_info.factor_group().size() << std::endl;
  for(auto const &op : supercell_sym_info.factor_group()) {
    log() << i << ": (" << op.index() << ") " << brief_description(op, supercell_sym_info.supercell_lattice()) << std::endl;
    ++i;
  }
  i = 0;
  for(auto const &cluster : scel_generators) {
    log() << "i: " << i << std::endl;
    printer.print(cluster, log());
    ++i;
  }

  /// Orbits based on primitive configuration factor group
  ScelPeriodicSymCompare<IntegralCluster> scel_sym_compare {
    shared_prim,
    supercell_sym_info.transformation_matrix_to_super(),
    shared_prim->lattice().tol()};
  auto prim_config_fg = configuration_L12.factor_group();
  std::vector<Permutation> prim_config_inverse_permutations = make_inverse_permutations(
                                                                prim_config_fg.begin(),
                                                                prim_config_fg.end());
  std::vector<IntegralCluster> orbit_generators_prim_config;
  make_suborbit_generators(
    supercell_sym_info.permute_begin(),
    supercell_sym_info.permute_end(),
    prim_config_fg.begin(),
    prim_config_fg.end(),
    scel_generators.begin(),
    scel_generators.end(),
    scel_sym_compare,
    std::back_inserter(orbit_generators_prim_config));

  i = 0;
  log() << "prim_config_fg.size(): " << prim_config_fg.size() << std::endl;
  for(auto const &permute_it : prim_config_fg) {
    auto op = permute_it->sym_op();
    log() << i << ": (" << op.index() << ") " << brief_description(op, supercell_sym_info.supercell_lattice()) << std::endl;
    ++i;
  }
  i = 0;
  for(auto const &cluster : orbit_generators_prim_config) {
    log() << "i: " << i << std::endl;
    printer.print(cluster, log());
    auto cluster_site_indices = make_cluster_site_indices(cluster, supercell_sym_info);
    auto canonical_indices = make_canonical_cluster_site_indices(
                               prim_config_inverse_permutations.begin(),
                               prim_config_inverse_permutations.end(),
                               cluster_site_indices);
    log() << "sites: " << jsonParser(cluster_site_indices)
          << "  canonical: " << jsonParser(canonical_indices) << std::endl;
    ++i;
  }

  /// Orbits based on non-primitive configuration factor group, under periodic boundary conditions
  auto superconfig_fg = configuration_L12_big.factor_group();
  auto supersupercell_sym_info = configuration_L12_big.supercell().sym_info();
  std::vector<Permutation> super_config_inverse_permutations = make_inverse_permutations(
                                                                 superconfig_fg.begin(),
                                                                 superconfig_fg.end());

  auto orbit_generators_superconfig_pbc = make_orbit_generators_under_periodic_boundary_conditions(
                                            supersupercell_sym_info,
                                            super_config_inverse_permutations.begin(),
                                            super_config_inverse_permutations.end(),
                                            orbit_generators_prim_config.begin(),
                                            orbit_generators_prim_config.end());

  i = 0;
  for(auto const &cluster_site_indices : orbit_generators_superconfig_pbc) {
    log() << "i: " << i << "  sites: " << jsonParser(cluster_site_indices) << std::endl;
    ++i;
  }


  // // Create ConfigEnumInput that are copies of `reference_config_enum_input` with orbit prototype
  // // sites selected, for all orbits generated by "cluster_specs"
  // std::vector<ConfigEnumInput> with_cluster_sites = select_cluster_sites(
  //                                                     reference_config_enum_input,
  //                                                     cluster_specs);
  // EXPECT_EQ(with_cluster_sites.size(), 5); // Expect 5 results
  // EXPECT_EQ(count_selected(with_cluster_sites, {0, 0, 0}), 1); // expect 1 unselected
  // EXPECT_EQ(count_selected(with_cluster_sites, {1, 0, 0}), 1); // expect 1 with 1 "A" site selected
  // EXPECT_EQ(count_selected(with_cluster_sites, {0, 1, 0}), 1); // expect 1 with 1 "B" site selected
  // EXPECT_EQ(count_selected(with_cluster_sites, {1, 1, 0}), 1); // expect 1 with 1 "A" and 1 "B" site selected
  // EXPECT_EQ(count_selected(with_cluster_sites, {2, 0, 0}), 1); // expect 1 with 2 "A" site selected

}

// TEST_F(ClusterSitesSelectorTest, Test2) {
//
//   // Example, generating Configuration with unique cluster sites selected:
//   // (starting from "within_scel" clusters)
//   // - L12 (A3B1) in conventional FCC supercell (4 atoms)
//
//   ConfigEnumInput reference_config_enum_input {configuration_L12_3x3x3, {}};
//
//   double a = configuration_L12.supercell().lattice()[0].norm();
//   std::vector<double> max_length {
//     0, // null cluster
//     0, // point clusters
//     3 * a * sqrt(3.) // a * sqrt(2.0) / 2.0 + TOL // 1NN pair clusters
//   };
//
//   Supercell const &supercell = reference_config_enum_input.configuration().supercell();
//   WithinScelMaxLengthClusterSpecs cluster_specs {
//     shared_prim,
//     &supercell.sym_info(),
//     make_generating_group(reference_config_enum_input),
//     alloy_sites_filter,
//     max_length};
//
//   // Create ConfigEnumInput that are copies of `reference_config_enum_input` with orbit prototype
//   // sites selected, for all orbits generated by "cluster_specs"
//   std::vector<ConfigEnumInput> with_cluster_sites = select_cluster_sites(
//                                                       reference_config_enum_input,
//                                                       cluster_specs);
//
//   EXPECT_EQ(with_cluster_sites.size(), 5); // Expect 3 results
//   EXPECT_EQ(count_selected(with_cluster_sites, {0, 0, 0}), 1); // expect 1 unselected
//   EXPECT_EQ(count_selected(with_cluster_sites, {1, 0, 0}), 1); // expect 1 with 1 "A" site selected
//   EXPECT_EQ(count_selected(with_cluster_sites, {0, 1, 0}), 1); // expect 1 with 1 "B" site selected
//   EXPECT_EQ(count_selected(with_cluster_sites, {1, 1, 0}), 1); // expect 1 with 1 "A" and 1 "B" site selected
//   EXPECT_EQ(count_selected(with_cluster_sites, {2, 0, 0}), 1); // expect 1 with 2 "A" site selected
//
// }

namespace {

  /// Make transformation_matrix_to_super for conventional 4 atom FCC cell, from the 1 atom primitive
  Eigen::Matrix3i _fcc_conventional_transf_mat() {
    Eigen::Matrix3i transf_mat;
    transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    return transf_mat;
  }

  /// Make L12 (A3B1), in conventional 4 atom FCC cell
  Configuration _make_configuration_L12(std::shared_ptr<CASM::Structure const> shared_prim) {
    auto conventional_fcc = std::make_shared<Supercell>(shared_prim, _fcc_conventional_transf_mat());
    Configuration configuration {conventional_fcc};
    configuration.configdof().occ(0) = 1;
    return configuration;
  }

  /// 3x3x3 super configuration of configuration_L12
  Configuration _make_configuration_L12_3x3x3(Configuration const &_configuration_L12) {
    auto conventional_fcc_3x3x3 = std::make_shared<Supercell>(
                                    _configuration_L12.supercell().shared_prim(),
                                    3 * _fcc_conventional_transf_mat());
    return fill_supercell(_configuration_L12, conventional_fcc_3x3x3);
  }

  /// NxNxN super configuration of configuration_L12
  Configuration _make_configuration_L12_big(long N, Configuration const &_configuration_L12) {
    auto conventional_fcc_big = std::make_shared<Supercell>(
                                  _configuration_L12.supercell().shared_prim(),
                                  N * _fcc_conventional_transf_mat());
    return fill_supercell(_configuration_L12, conventional_fcc_big);
  }

  /// Check how many ConfigEnumInput have a matching number of each type of occupant on the selected sites
  ///
  /// Example, FCC ternary prim:
  /// - Count how many ConfigEnumInput in `with_cluster_sites` have 1 site occupied by an A atom,
  ///   2 sites occupied by a B atom, and 0 sites occupied by a C atom selected:
  ///       `count_selected(with_cluster_sites, {1, 2, 0})`
  Index count_selected(std::vector<ConfigEnumInput> const &with_cluster_sites, std::vector<Index> expected_occ_count) {
    Index matching = 0;
    for(ConfigEnumInput const &config_with_cluster : with_cluster_sites) {
      Configuration const &config = config_with_cluster.configuration();
      std::vector<Index> occ_count(expected_occ_count.size(), 0);
      for(Index site_index : config_with_cluster.sites()) {
        occ_count[config.occ(site_index)]++;
      }
      if(occ_count == expected_occ_count) {
        matching++;
      }
    }
    return matching;
  }

}
