#include "autotools.hh"
#include "casm/app/AppIO_impl.hh"  // for orbit printing
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
#include "gtest/gtest.h"

using namespace CASM;
using namespace test;

// helper functions
namespace {

/// Make transformation_matrix_to_super for conventional 4 atom FCC cell, from
/// the 1 atom primitive
Eigen::Matrix3l _fcc_conventional_transf_mat();

/// Make L12 (A3B1), in conventional 4 atom FCC cell
Configuration _make_configuration_L12(
    std::shared_ptr<CASM::Structure const> shared_prim);

/// NxNxN super configuration of configuration_L12
Configuration _make_super_configuration_L12(
    long N, Configuration const &_configuration_L12);

/// Check how many clusters have a matching number of each type of occupant on
/// the selected sites
Index count_selected(
    Configuration const &configuration,
    std::set<std::set<Index>, ClusterSiteIndicesCompare> const &generators,
    std::vector<Index> expected_occ_count);
}  // namespace

class ClusterSitesSelectorTest : public testing::Test {
 protected:
  std::shared_ptr<CASM::Structure const> shared_prim;
  Configuration
      configuration_L12;  // L12 (A3B1), in conventional 4 atom FCC cell
  Configuration
      superconfiguration_L12;  // superconfiguration of configuration_L12

  ClusterSitesSelectorTest()
      : shared_prim(
            std::make_shared<CASM::Structure const>(test::FCC_ternary_prim())),
        configuration_L12(_make_configuration_L12(shared_prim)),
        superconfiguration_L12(
            _make_super_configuration_L12(5, configuration_L12)) {
    EXPECT_EQ(configuration_L12.size(), 4);
    EXPECT_EQ(superconfiguration_L12.size(), pow(5, 3) * 4);
  }
};

TEST_F(ClusterSitesSelectorTest, Test1) {
  ScopedNullLogging logging;

  // Example, generating Configuration with unique cluster sites selected:
  // (starting from "periodic_max_length" clusters)
  // - L12 (A3B1) in conventional FCC supercell (4 atoms)

  ConfigEnumInput reference_config_enum_input{configuration_L12, {}};

  // double a = shared_prim->lattice()[0].norm();
  double a = configuration_L12.supercell().lattice()[0].norm();

  std::vector<double> max_length{
      0,            // null cluster
      0,            // point clusters
      2 * a + TOL,  // a * sqrt(2.0) / 2.0 + TOL // 1NN pair clusters
      2 * a + TOL   // a * sqrt(2.0) / 2.0 + TOL // 1NN pair clusters
  };

  PeriodicMaxLengthClusterSpecs cluster_specs{
      shared_prim, shared_prim->factor_group(), alloy_sites_filter,
      max_length};  // max length, for point clusters only (null orbit branch,
                    // point orbit branch)
  auto orbits = cluster_specs.make_periodic_orbits(log());

  SupercellSymInfo const &supercell_sym_info =
      configuration_L12.supercell().sym_info();
  std::vector<IntegralCluster> scel_generators;
  make_suborbit_generators(
      shared_prim->factor_group().begin(), shared_prim->factor_group().end(),
      supercell_sym_info.factor_group().begin(),
      supercell_sym_info.factor_group().end(),
      prototype_iterator(orbits.begin()), prototype_iterator(orbits.end()),
      cluster_specs.sym_compare, std::back_inserter(scel_generators));

  // // Uncomment to print scel orbit generators
  // {
  //   OrbitPrinterOptions printer_options;
  //   printer_options.coord_type = CART;
  //   printer_options.delim = char();
  //   printer_options.orbit_print_mode = ORBIT_PRINT_MODE::PROTO;
  //   Printer<IntegralCluster> printer {printer_options};
  //
  //   Index i;
  //   i = 0;
  //   log() << "supercell_sym_info.factor_group().size(): "
  //         << supercell_sym_info.factor_group().size() << std::endl;
  //   for(auto const &op : supercell_sym_info.factor_group()) {
  //     log() << i << ": (" << op.index() << ") " << brief_description(op,
  //     supercell_sym_info.supercell_lattice()) << std::endl;
  //     ++i;
  //   }
  //   i = 0;
  //   for(auto const &cluster : scel_generators) {
  //     log() << "i: " << i << std::endl;
  //     printer.print(cluster, log());
  //     ++i;
  //   }
  // }

  /// Orbits based on primitive configuration factor group
  ScelPeriodicSymCompare<IntegralCluster> scel_sym_compare{
      shared_prim, supercell_sym_info.transformation_matrix_to_super(),
      shared_prim->lattice().tol()};
  auto prim_config_fg = configuration_L12.factor_group();
  std::vector<Permutation> prim_config_inverse_permutations =
      make_inverse_permutations(prim_config_fg.begin(), prim_config_fg.end());
  std::vector<IntegralCluster> orbit_generators_prim_config;
  make_suborbit_generators(
      supercell_sym_info.permute_begin(), supercell_sym_info.permute_end(),
      prim_config_fg.begin(), prim_config_fg.end(), scel_generators.begin(),
      scel_generators.end(), scel_sym_compare,
      std::back_inserter(orbit_generators_prim_config));

  // // Uncomment to print primitive L12 orbit generators
  // {
  //   OrbitPrinterOptions printer_options;
  //   printer_options.coord_type = CART;
  //   printer_options.delim = char();
  //   printer_options.orbit_print_mode = ORBIT_PRINT_MODE::PROTO;
  //   Printer<IntegralCluster> printer {printer_options};
  //
  //   Index i;
  //   i = 0;
  //   log() << "prim_config_fg.size(): " << prim_config_fg.size() << std::endl;
  //   for(auto const &permute_it : prim_config_fg) {
  //     auto op = permute_it->sym_op();
  //     log() << i << ": (" << op.index() << ") " << brief_description(op,
  //     supercell_sym_info.supercell_lattice()) << std::endl;
  //     ++i;
  //   }
  //   i = 0;
  //   for(auto const &cluster : orbit_generators_prim_config) {
  //     log() << "i: " << i << std::endl;
  //     printer.print(cluster, log());
  //     auto cluster_site_indices = make_cluster_site_indices(cluster,
  //     supercell_sym_info); auto canonical_indices =
  //     make_canonical_cluster_site_indices(
  //                                prim_config_inverse_permutations.begin(),
  //                                prim_config_inverse_permutations.end(),
  //                                cluster_site_indices);
  //     log() << "sites: " << jsonParser(cluster_site_indices)
  //           << "  canonical: " << jsonParser(canonical_indices) << std::endl;
  //     ++i;
  //   }
  // }

  /// Orbits based on primitive configuration factor group, under periodic
  /// boundary conditions
  auto unitconfig_fg = configuration_L12.factor_group();
  auto unitsupercell_sym_info = configuration_L12.supercell().sym_info();
  std::vector<Permutation> unit_config_inverse_permutations =
      make_inverse_permutations(unitconfig_fg.begin(), unitconfig_fg.end());

  std::set<std::set<Index>, ClusterSiteIndicesCompare>
      orbit_generators_unitconfig_pbc =
          make_orbit_generators_under_periodic_boundary_conditions(
              unitsupercell_sym_info, unit_config_inverse_permutations.begin(),
              unit_config_inverse_permutations.end(),
              orbit_generators_prim_config.begin(),
              orbit_generators_prim_config.end());

  // check number of clusters occupying {#A, #B, #C} sites:
  EXPECT_EQ(count_selected(configuration_L12, orbit_generators_unitconfig_pbc,
                           {0, 0, 0}),
            1);
  EXPECT_EQ(count_selected(configuration_L12, orbit_generators_unitconfig_pbc,
                           {1, 0, 0}),
            1);
  EXPECT_EQ(count_selected(configuration_L12, orbit_generators_unitconfig_pbc,
                           {0, 1, 0}),
            1);
  EXPECT_EQ(count_selected(configuration_L12, orbit_generators_unitconfig_pbc,
                           {1, 1, 0}),
            1);
  EXPECT_EQ(count_selected(configuration_L12, orbit_generators_unitconfig_pbc,
                           {2, 0, 0}),
            1);
  EXPECT_EQ(count_selected(configuration_L12, orbit_generators_unitconfig_pbc,
                           {2, 1, 0}),
            1);
  EXPECT_EQ(count_selected(configuration_L12, orbit_generators_unitconfig_pbc,
                           {3, 0, 0}),
            1);

  /// Orbits based on non-primitive configuration factor group, under periodic
  /// boundary conditions
  auto superconfig_fg = superconfiguration_L12.factor_group();
  auto supersupercell_sym_info = superconfiguration_L12.supercell().sym_info();
  std::vector<Permutation> super_config_inverse_permutations =
      make_inverse_permutations(superconfig_fg.begin(), superconfig_fg.end());

  std::set<std::set<Index>, ClusterSiteIndicesCompare>
      orbit_generators_superconfig_pbc =
          make_orbit_generators_under_periodic_boundary_conditions(
              supersupercell_sym_info,
              super_config_inverse_permutations.begin(),
              super_config_inverse_permutations.end(),
              orbit_generators_prim_config.begin(),
              orbit_generators_prim_config.end());

  // Uncomment to print superconfiguration L12 orbit generators (cluster site
  // indices as std::set<Index>)
  {
    Index i = 0;
    for (auto const &cluster_site_indices : orbit_generators_superconfig_pbc) {
      log() << "i: " << i << "  sites: " << jsonParser(cluster_site_indices)
            << std::endl;
      ++i;
    }
  }

  // when the supercell is large enough, pbc do not cause aliasing, and the
  // number of generators is unchanged
  EXPECT_EQ(orbit_generators_superconfig_pbc.size(),
            orbit_generators_prim_config.size());
}

namespace {

/// Make transformation_matrix_to_super for conventional 4 atom FCC cell, from
/// the 1 atom primitive
Eigen::Matrix3l _fcc_conventional_transf_mat() {
  Eigen::Matrix3l transf_mat;
  transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
  return transf_mat;
}

/// Make L12 (A3B1), in conventional 4 atom FCC cell
Configuration _make_configuration_L12(
    std::shared_ptr<CASM::Structure const> shared_prim) {
  auto conventional_fcc =
      std::make_shared<Supercell>(shared_prim, _fcc_conventional_transf_mat());
  Configuration configuration{conventional_fcc};
  configuration.configdof().occ(0) = 1;
  return configuration;
}

/// NxNxN super configuration of configuration_L12
Configuration _make_super_configuration_L12(
    long N, Configuration const &_configuration_L12) {
  auto conventional_fcc_supercell =
      std::make_shared<Supercell>(_configuration_L12.supercell().shared_prim(),
                                  N * _fcc_conventional_transf_mat());
  return fill_supercell(_configuration_L12, conventional_fcc_supercell);
}

/// Check how many ConfigEnumInput have a matching number of each type of
/// occupant on the selected sites
///
/// Example, FCC ternary prim:
/// - Count how many generators have 1 site occupied by an A atom, 2 sites
/// occupied by a B atom,
///   and 0 sites occupied by a C atom selected:
///       `count_selected(configuration, generators, {1, 2, 0})`
Index count_selected(
    Configuration const &configuration,
    std::set<std::set<Index>, ClusterSiteIndicesCompare> const &generators,
    std::vector<Index> expected_occ_count) {
  Index matching = 0;
  for (std::set<Index> const &cluster_site_indices : generators) {
    std::vector<Index> occ_count(expected_occ_count.size(), 0);
    for (Index site_index : cluster_site_indices) {
      occ_count[configuration.occ(site_index)]++;
    }
    if (occ_count == expected_occ_count) {
      matching++;
    }
  }
  return matching;
}

}  // namespace
