#include "autotools.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "casm/clusterography/io/OrbitPrinter_impl.hh"
#include "casm/clusterography/io/stream/IntegralCluster_stream_io.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/SimpleOrbit_impl.hh"
#include "casm/symmetry/SubOrbits_impl.hh"
#include "casm/symmetry/SupercellSymInfo.hh"
#include "crystallography/TestStructures.hh"
#include "gtest/gtest.h"

using namespace CASM;
using namespace test;

namespace {
Eigen::Matrix3l _fcc_conventional_transf_mat() {
  Eigen::Matrix3l transf_mat;
  transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
  return transf_mat;
}
}  // namespace

namespace test {

void expect_excluded_point_cluster(
    std::vector<Orbit<LocalSymCompare<IntegralCluster>>> const &orbits,
    xtal::UnitCellCoord const &unitcellcoord);

void expect_included_point_cluster(
    std::vector<Orbit<LocalSymCompare<IntegralCluster>>> const &orbits,
    xtal::UnitCellCoord const &unitcellcoord, int multiplicity);
}  // namespace test

class ClusterSpecsTest : public testing::Test {
 protected:
  std::shared_ptr<CASM::Structure const> shared_prim;

  ClusterSpecsTest()
      : shared_prim(
            std::make_shared<CASM::Structure const>(test::FCC_ternary_prim())) {
  }
};

TEST_F(ClusterSpecsTest, PeriodicTest) {
  // primitive FCC clusters

  // a = primitive FCC unit cell lattice vector length
  double a = shared_prim->lattice()[0].norm();

  std::vector<double> max_length{
      0,       // null cluster
      0,       // point clusters
      a + TOL  // 1NN pair clusters only
  };

  PeriodicMaxLengthClusterSpecs cluster_specs{
      shared_prim, shared_prim->factor_group(), alloy_sites_filter, max_length};

  auto orbits = cluster_specs.make_periodic_orbits(null_log());

  EXPECT_EQ(orbits.size(), 3);
  EXPECT_EQ(orbits[0].prototype().size(), 0);  // 1 null cluster
  EXPECT_EQ(orbits[0].size(), 1);
  EXPECT_EQ(orbits[1].prototype().size(), 1);  // 1 point cluster
  EXPECT_EQ(orbits[1].size(), 1);
  EXPECT_EQ(orbits[2].prototype().size(),
            2);  // 6 1NN pair clusters (per site, no double counting)
  EXPECT_EQ(orbits[2].size(), 6);
}

TEST_F(ClusterSpecsTest, ScelPeriodicTest) {
  // Make "scel_periodic" orbit generators from "prim_periodic" orbits

  // 1) make periodic orbits: (1NN pair clusters only)

  // a = primitive FCC unit cell lattice vector length
  double a = shared_prim->lattice()[0].norm();

  std::vector<double> max_length{
      0,       // null cluster
      0,       // point clusters
      a + TOL  // 1NN pair clusters only
  };

  PeriodicMaxLengthClusterSpecs cluster_specs{
      shared_prim, shared_prim->factor_group(), alloy_sites_filter, max_length};

  auto prim_periodic_orbits = cluster_specs.make_periodic_orbits(null_log());

  // // Uncomment to print the prim factor group:
  // {
  //   Index i = 0;
  //   log() << "shared_prim->factor_group().size(): "
  //         << shared_prim->factor_group().size() << std::endl;
  //   for(auto const &op : shared_prim->factor_group()) {
  //     log() << i++ << ": (" << op.index() << ") " << brief_description(op,
  //     shared_prim->lattice()) << std::endl;
  //   }
  // }

  // // Uncomment the following to print the prim_periodic_orbits:
  // // Construct OrbitPrinterOptions:
  // {
  //   OrbitPrinterOptions printer_options;
  //   printer_options.coord_type = CART;
  //   printer_options.delim = char();
  //   printer_options.orbit_print_mode = ORBIT_PRINT_MODE::PROTO;
  //   Printer<IntegralCluster> printer {printer_options};
  //
  //   print_clust(
  //     prim_periodic_orbits.begin(),
  //     prim_periodic_orbits.end(),
  //     log(),
  //     printer_options);
  // }

  // null, point, 1NN pair orbits
  EXPECT_EQ(prim_periodic_orbits.size(), 3);
  EXPECT_EQ(prim_periodic_orbits[0].prototype().size(), 0);
  EXPECT_EQ(prim_periodic_orbits[1].prototype().size(), 1);
  EXPECT_EQ(prim_periodic_orbits[2].prototype().size(), 2);
  EXPECT_EQ(prim_periodic_orbits[2].size(), 6);

  // 2) make a supercell (3x2x1 supercell of the conventional FCC unit cell)

  // conventional FCC (4 atom unit cell) clusters "within_scel"
  Eigen::Matrix3l T = _fcc_conventional_transf_mat();

  // U makes supercells of the conventional FCC unit cell
  Eigen::Matrix3l U;
  U << 3, 0, 0, 0, 2, 0, 0, 0, 1;

  // combine to construct prim -> supercell transformation matrix
  T = _fcc_conventional_transf_mat() * U;

  Lattice super_lattice = make_superlattice(shared_prim->lattice(), T);
  SupercellSymInfo supercell_sym_info =
      make_supercell_sym_info(*shared_prim, super_lattice);

  // 3) Find sub orbits due to prim factor group -> supercell factor group
  // symmetry breaking

  std::vector<IntegralCluster> suborbit_generators;
  make_suborbit_generators(
      shared_prim->factor_group().begin(), shared_prim->factor_group().end(),
      supercell_sym_info.factor_group().begin(),
      supercell_sym_info.factor_group().end(),
      prototype_iterator(prim_periodic_orbits.begin()),
      prototype_iterator(prim_periodic_orbits.end()), cluster_specs.sym_compare,
      std::back_inserter(suborbit_generators));

  // // Uncomment to print the supercell factor group:
  // {
  //   Index i = 0;
  //   log() << "supercell_sym_info.factor_group().size(): "
  //         << supercell_sym_info.factor_group().size() << std::endl;
  //   for(auto const &op : supercell_sym_info.factor_group()) {
  //     log() << i++ << ": (" << op.index() << ") " << brief_description(op,
  //     shared_prim->lattice()) << std::endl;
  //   }
  // }

  // // Uncomment the following to print the sub-orbit generators:
  // {
  //   OrbitPrinterOptions printer_options;
  //   printer_options.coord_type = CART;
  //   printer_options.delim = char();
  //   printer_options.orbit_print_mode = ORBIT_PRINT_MODE::PROTO;
  //   Printer<IntegralCluster> printer {printer_options};
  //
  //   Index i = 0;
  //   for(auto const &cluster : suborbit_generators) {
  //     log() << "i: " << i << std::endl;
  //     printer.print(cluster, log());
  //     ++i;
  //   }
  // }

  // 4) Generate "scel_periodic" orbits:
  std::vector<SimpleOrbit<CASM::ScelPeriodicSymCompare<IntegralCluster>>>
      scel_periodic_orbits;
  CASM::ScelPeriodicSymCompare<IntegralCluster> scel_periodic_sym_compare{
      shared_prim, supercell_sym_info.transformation_matrix_to_super(),
      shared_prim->lattice().tol()};
  for (auto const &generating_element : suborbit_generators) {
    scel_periodic_orbits.emplace_back(
        generating_element, supercell_sym_info.permute_begin(),
        supercell_sym_info.permute_end(), scel_periodic_sym_compare);
  }

  // // Uncomment the following to print the scel_periodic_orbits:
  // {
  //   OrbitPrinterOptions printer_options;
  //   printer_options.coord_type = CART;
  //   printer_options.delim = char();
  //   printer_options.orbit_print_mode = ORBIT_PRINT_MODE::PROTO;
  //   Printer<IntegralCluster> printer {printer_options};
  //
  //   print_clust(
  //     scel_periodic_orbits.begin(),
  //     scel_periodic_orbits.end(),
  //     log(),
  //     printer_options);
  // }

  // null, point, 1NN pair orbits
  EXPECT_EQ(scel_periodic_orbits.size(), 5);
  EXPECT_EQ(scel_periodic_orbits[0].prototype().size(), 0);
  EXPECT_EQ(scel_periodic_orbits[1].prototype().size(), 1);
  EXPECT_EQ(scel_periodic_orbits[2].prototype().size(), 2);
  EXPECT_EQ(scel_periodic_orbits[2].size(),
            48);  // 2 elements per unit cell * 24 unit cells
  EXPECT_EQ(scel_periodic_orbits[3].prototype().size(), 2);
  EXPECT_EQ(scel_periodic_orbits[3].size(),
            48);  // 2 elements per unit cell * 24 unit cells
  EXPECT_EQ(scel_periodic_orbits[4].prototype().size(), 2);
  EXPECT_EQ(scel_periodic_orbits[4].size(),
            48);  // 2 elements per unit cell * 24 unit cells
}

TEST_F(ClusterSpecsTest, LocalTest) {
  // primitive FCC: local point clusters about the 1NN pair cluster

  // a = primitive FCC unit cell lattice vector length
  double a = shared_prim->lattice()[0].norm();

  // 1NN phenomenal cluster: 1NN pair cluster: {0, 0, 0} & {a, a, 0}
  IntegralCluster phenomenal_cluster{*shared_prim};
  phenomenal_cluster.elements().emplace_back(0, xtal::UnitCell{0, 0, 0});
  phenomenal_cluster.elements().emplace_back(0, xtal::UnitCell{0, 0, 1});

  // Make the group that leaves the phenomenal cluster invariant:
  //
  // Note: Use `PrimPeriodicSymCompare` here to include all factor group
  // operation + translation
  //       combinations in the cluster's invariant group. The
  //       `spatial_transform` necessary for equivalence is included in the
  //       invariant subgroup operation.
  PrimPeriodicSymCompare<IntegralCluster> prim_periodic_sym_compare{
      shared_prim, shared_prim->lattice().tol()};
  auto phenomenal_cluster_group =
      make_invariant_subgroup(phenomenal_cluster, shared_prim->factor_group(),
                              prim_periodic_sym_compare);

  // The value max_length[b], is the max site-to-site distance within a cluster
  // for that cluster to be included in branch b. The b==0 value is ignored.
  std::vector<double> max_length{
      0,  // null cluster
      0   // point clusters
  };

  // For a site to be added to clusters in branch b it must be a distance less
  // than cutoff_radius[b] to any site in the phenomenal cluster. The b==0 value
  // is ignored.
  std::vector<double> cutoff_radius{
      0,       // null cluster
      a + TOL  // point clusters
  };

  // If true, local clusters include phenomenal_cluster sites; otherwise they do
  // not
  bool include_phenomenal_sites = false;

  LocalMaxLengthClusterSpecs cluster_specs{shared_prim,
                                           phenomenal_cluster_group,
                                           phenomenal_cluster,
                                           alloy_sites_filter,
                                           max_length,
                                           cutoff_radius,
                                           include_phenomenal_sites};

  auto orbits = cluster_specs.make_local_orbits(null_log());

  // // Uncomment the following to print the local orbits:
  // {
  //   OrbitPrinterOptions printer_options;
  //   printer_options.coord_type = CART;
  //   printer_options.sym_info_opt = SymInfoOptions
  //   {printer_options.coord_type}; printer_options.delim = char();
  //   printer_options.orbit_print_mode = ORBIT_PRINT_MODE::PROTO;
  //   Printer<IntegralCluster> printer {printer_options};
  //
  //   COORD_MODE mode {CART};
  //
  //   log() << "phenomenal:\n";
  //   printer.print(
  //     phenomenal_cluster,
  //     log());
  //   log() << std::endl;
  //
  //   log() << "factor group:\n";
  //   brief_description(log(), shared_prim->factor_group(),
  //   shared_prim->lattice(), printer_options.sym_info_opt); log() <<
  //   std::endl;
  //
  //   log() << "invariant group:\n";
  //   brief_description(log(), phenomenal_cluster_group,
  //   shared_prim->lattice(), printer_options.sym_info_opt); log() <<
  //   std::endl;
  //
  //   print_clust(
  //     orbits.begin(),
  //     orbits.end(),
  //     log(),
  //     printer_options);
  // }

  EXPECT_EQ(orbits.size(), 5);
  EXPECT_EQ(orbits[0].prototype().size(), 0);  // 1 null cluster
  EXPECT_EQ(orbits[0].size(), 1);

  // expect point cluster {b, i, j, k} to not be found in any orbit
  expect_excluded_point_cluster(orbits, {0, 0, 0, 0});

  // expect point cluster {b, i, j, k} to be found in an orbit with specified
  // multiplicity
  expect_included_point_cluster(orbits, {0, 1, 0, -1}, 8);
  expect_included_point_cluster(orbits, {0, -1, 1, 1}, 4);
  expect_included_point_cluster(orbits, {0, 0, 1, 0}, 4);
  expect_included_point_cluster(orbits, {0, 0, 0, -1}, 2);
}

namespace test {

// Check that point cluster consisiting of `unitcellcoord` is not found in
// orbits
void expect_excluded_point_cluster(
    std::vector<Orbit<LocalSymCompare<IntegralCluster>>> const &orbits,
    xtal::UnitCellCoord const &unitcellcoord) {
  IntegralCluster test{orbits.begin()->prototype().prim()};
  test.elements().push_back(unitcellcoord);
  for (auto const &orbit : orbits) {
    EXPECT_EQ(orbit.contains(test), false);
  }
}

// Check that point cluster consisiting of `unitcellcoord` is found in orbits
// with `multiplicity`
void expect_included_point_cluster(
    std::vector<Orbit<LocalSymCompare<IntegralCluster>>> const &orbits,
    xtal::UnitCellCoord const &unitcellcoord, int multiplicity) {
  IntegralCluster test{orbits.begin()->prototype().prim()};
  test.elements().push_back(unitcellcoord);
  for (auto const &orbit : orbits) {
    if (orbit.contains(test)) {
      EXPECT_EQ(orbit.size(), multiplicity);
      return;
    }
  }
  EXPECT_EQ(true, false) << "point cluster expected to be found was not found";
}
}  // namespace test
