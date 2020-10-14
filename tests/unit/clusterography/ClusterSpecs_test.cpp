#include "gtest/gtest.h"
#include "autotools.hh"
#include "casm/app/AppIO_impl.hh" // for orbit printing
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
// #include "casm/clusterography/SupercellClusterOrbits_impl.hh"
#include "casm/clusterography/io/stream/IntegralCluster_stream_io.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/SimpleOrbit_impl.hh"
#include "casm/symmetry/SubOrbits_impl.hh"
#include "casm/symmetry/SupercellSymInfo_impl.hh"
#include "crystallography/TestStructures.hh"

using namespace CASM;
using namespace test;

namespace {
  Eigen::Matrix3l _fcc_conventional_transf_mat() {
    Eigen::Matrix3l transf_mat;
    transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    return transf_mat;
  }
}

class ClusterSpecsTest : public testing::Test {
protected:

  std::shared_ptr<CASM::Structure const> shared_prim;

  ClusterSpecsTest():
    shared_prim(std::make_shared<CASM::Structure const>(test::FCC_ternary_prim())) {}

};

TEST_F(ClusterSpecsTest, PeriodicTest) {

  // primitive FCC clusters

  // a = primitive FCC unit cell lattice vector length
  double a = shared_prim->lattice()[0].norm();

  std::vector<double> max_length {
    0, // null cluster
    0, // point clusters
    a + TOL // 1NN pair clusters only
  };

  PeriodicMaxLengthClusterSpecs cluster_specs {
    shared_prim,
    shared_prim->factor_group(),
    alloy_sites_filter,
    max_length};

  auto orbits = cluster_specs.make_periodic_orbits(null_log());

  EXPECT_EQ(orbits.size(), 3);
  EXPECT_EQ(orbits[0].prototype().size(), 0); // 1 null cluster
  EXPECT_EQ(orbits[0].size(), 1);
  EXPECT_EQ(orbits[1].prototype().size(), 1); // 1 point cluster
  EXPECT_EQ(orbits[1].size(), 1);
  EXPECT_EQ(orbits[2].prototype().size(), 2); // 6 1NN pair clusters (per site, no double counting)
  EXPECT_EQ(orbits[2].size(), 6);
}

TEST_F(ClusterSpecsTest, ScelPeriodicTest) {

  // Make "scel_periodic" orbit generators from "prim_periodic" orbits

  // 1) make periodic orbits: (1NN pair clusters only)

  // a = primitive FCC unit cell lattice vector length
  double a = shared_prim->lattice()[0].norm();

  std::vector<double> max_length {
    0, // null cluster
    0, // point clusters
    a + TOL // 1NN pair clusters only
  };

  PeriodicMaxLengthClusterSpecs cluster_specs {
    shared_prim,
    shared_prim->factor_group(),
    alloy_sites_filter,
    max_length};

  auto prim_periodic_orbits = cluster_specs.make_periodic_orbits(null_log());

  // Uncomment to print the prim factor group:
  {
    Index i = 0;
    log() << "shared_prim->factor_group().size(): "
          << shared_prim->factor_group().size() << std::endl;
    for(auto const &op : shared_prim->factor_group()) {
      log() << i++ << ": (" << op.index() << ") " << brief_description(op, shared_prim->lattice()) << std::endl;
    }
  }

  // Uncomment the following to print the prim_periodic_orbits:
  // Construct OrbitPrinterOptions:
  {
    OrbitPrinterOptions printer_options;
    printer_options.coord_type = CART;
    printer_options.delim = char();
    printer_options.orbit_print_mode = ORBIT_PRINT_MODE::FULL;
    Printer<IntegralCluster> printer {printer_options};

    print_clust(
      prim_periodic_orbits.begin(),
      prim_periodic_orbits.end(),
      log(),
      printer_options);
  }

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
  U << 3, 0, 0,
  0, 2, 0,
  0, 0, 1;

  // combine to construct prim -> supercell transformation matrix
  T = _fcc_conventional_transf_mat() * U;

  Lattice super_lattice = make_superlattice(shared_prim->lattice(), T);
  SupercellSymInfo supercell_sym_info = make_supercell_sym_info(*shared_prim, super_lattice);

  // 3) Find sub orbits due to prim factor group -> supercell factor group symmetry breaking

  std::vector<IntegralCluster> suborbit_generators;
  make_suborbit_generators(
    shared_prim->factor_group().begin(),
    shared_prim->factor_group().end(),
    supercell_sym_info.factor_group().begin(),
    supercell_sym_info.factor_group().end(),
    prototype_iterator(prim_periodic_orbits.begin()),
    prototype_iterator(prim_periodic_orbits.end()),
    cluster_specs.sym_compare,
    std::back_inserter(suborbit_generators));

  // Uncomment to print the supercell factor group:
  {
    Index i = 0;
    log() << "supercell_sym_info.factor_group().size(): "
          << supercell_sym_info.factor_group().size() << std::endl;
    for(auto const &op : supercell_sym_info.factor_group()) {
      log() << i++ << ": (" << op.index() << ") " << brief_description(op, shared_prim->lattice()) << std::endl;
    }
  }

  // Uncomment the following to print the sub-orbit generators:
  {
    OrbitPrinterOptions printer_options;
    printer_options.coord_type = CART;
    printer_options.delim = char();
    printer_options.orbit_print_mode = ORBIT_PRINT_MODE::PROTO;
    Printer<IntegralCluster> printer {printer_options};

    Index i = 0;
    for(auto const &cluster : suborbit_generators) {
      log() << "i: " << i << std::endl;
      printer.print(cluster, log());
      ++i;
    }
  }

  // 4) Generate "scel_periodic" orbits:

  std::vector<SimpleOrbit<CASM::ScelPeriodicSymCompare<IntegralCluster>>> scel_periodic_orbits;
  CASM::ScelPeriodicSymCompare<IntegralCluster> scel_periodic_sym_compare {
    shared_prim,
    supercell_sym_info.transformation_matrix_to_super(),
    shared_prim->lattice().tol()};
  for(auto const &generating_element : suborbit_generators) {
    scel_periodic_orbits.emplace_back(generating_element,
                                      supercell_sym_info.permute_begin(),
                                      supercell_sym_info.permute_end(),
                                      scel_periodic_sym_compare);
  }

  // Uncomment the following to print the scel_periodic_orbits:
  {
    OrbitPrinterOptions printer_options;
    printer_options.coord_type = CART;
    printer_options.delim = char();
    printer_options.orbit_print_mode = ORBIT_PRINT_MODE::FULL;
    Printer<IntegralCluster> printer {printer_options};

    print_clust(
      scel_periodic_orbits.begin(),
      scel_periodic_orbits.end(),
      log(),
      printer_options);
  }

  // null, point, 1NN pair orbits
  EXPECT_EQ(scel_periodic_orbits.size(), 5);
  EXPECT_EQ(scel_periodic_orbits[0].prototype().size(), 0);
  EXPECT_EQ(scel_periodic_orbits[1].prototype().size(), 1);
  EXPECT_EQ(scel_periodic_orbits[2].prototype().size(), 2);
  EXPECT_EQ(scel_periodic_orbits[2].size(), 6);
  EXPECT_EQ(scel_periodic_orbits[3].prototype().size(), 2);
  EXPECT_EQ(scel_periodic_orbits[3].size(), 6);
  EXPECT_EQ(scel_periodic_orbits[4].prototype().size(), 2);
  EXPECT_EQ(scel_periodic_orbits[4].size(), 6);
}

// TEST_F(ClusterSpecsTest, WithinScelTest) {
//
//   // conventional FCC (4 atom unit cell) clusters "within_scel"
//   Eigen::Matrix3l T = _fcc_conventional_transf_mat();
//   Lattice super_lattice = make_superlattice(shared_prim->lattice(), T);
//   SupercellSymInfo supercell_sym_info = make_supercell_sym_info(*shared_prim, super_lattice);
//
//   // a = primitive FCC unit cell lattice vector length
//   double a = shared_prim->lattice()[0].norm();
//
//   std::vector<double> max_length {
//     0, // null cluster
//     0, // point clusters
//     a + TOL // 1NN pair clusters only
//   };
//
//   WithinScelMaxLengthClusterSpecs cluster_specs {
//     shared_prim,
//     &supercell_sym_info,
//     make_generating_group(supercell_sym_info),
//     alloy_sites_filter,
//     max_length};
//
//   // auto orbits = cluster_specs.make_within_scel_orbits(log());
//   //
//   // // Construct OrbitPrinterOptions:
//   // OrbitPrinterOptions printer_options;
//   // printer_options.coord_type = INTEGRAL;
//   // // printer_options.print_equivalence_map = true;
//   // // printer_options.print_invariant_group = true;
//   //
//   // // Print prototype clusters and all equivalent clusters
//   // printer_options.orbit_print_mode = ORBIT_PRINT_MODE::PROTO;
//   // print_clust(orbits.begin(), orbits.end(), log(), printer_options);
//   //
//   // EXPECT_EQ(orbits.size(), 3);
//   //
//   // EXPECT_EQ(orbits[0].prototype().size(), 0); // 1 null cluster
//   // EXPECT_EQ(orbits[0].size(), 1);
//   // EXPECT_EQ(orbits[1].prototype().size(), 1); // 4 point clusters
//   // EXPECT_EQ(orbits[1].size(), 4);
//   // EXPECT_EQ(orbits[2].prototype().size(), 2); // 6 1NN pair clusters ("within" conventional unit cell, no double counting)
//   // EXPECT_EQ(orbits[2].size(), 6);
//
//   auto orbit_generators = cluster_specs.make_within_scel_orbit_generators(log());
//
//   jsonParser json;
//   json["orbit_generators"].put_array();
//   for(auto const &cluster_site_indices: orbit_generators) {
//     json["orbit_generators"].push_back(cluster_site_indices);
//   }
//   std::cout << json << std::endl;
//
// }
//
// TEST_F(ClusterSpecsTest, WithinScelTest2) {
//
//   // 2x2x2 conventional FCC (4 atom unit cell / 32 atom supercell) clusters "within_scel"
//   Eigen::Matrix3l T = 2 * _fcc_conventional_transf_mat();
//   Lattice super_lattice = make_superlattice(shared_prim->lattice(), T);
//   SupercellSymInfo supercell_sym_info = make_supercell_sym_info(*shared_prim, super_lattice);
//
//   // a = primitive FCC unit cell lattice vector length
//   double a = shared_prim->lattice()[0].norm();
//
//   std::vector<double> max_length {
//     0, // null cluster
//     0, // point clusters
//     6*a + TOL // 1NN pair clusters only
//   };
//
//   WithinScelMaxLengthClusterSpecs cluster_specs {
//     shared_prim,
//     &supercell_sym_info,
//     make_generating_group(supercell_sym_info),
//     alloy_sites_filter,
//     max_length};
//
//   // auto orbits = cluster_specs.make_within_scel_orbits(log());
//   //
//   // // Construct OrbitPrinterOptions:
//   // OrbitPrinterOptions printer_options;
//   // printer_options.coord_type = INTEGRAL;
//   // // printer_options.print_equivalence_map = true;
//   // // printer_options.print_invariant_group = true;
//   //
//   // // Print prototype clusters and all equivalent clusters
//   // printer_options.orbit_print_mode = ORBIT_PRINT_MODE::PROTO;
//   // print_clust(orbits.begin(), orbits.end(), log(), printer_options);
//   //
//   // EXPECT_EQ(orbits.size(), 3);
//   //
//   // EXPECT_EQ(orbits[0].prototype().size(), 0); // 1 null cluster
//   // EXPECT_EQ(orbits[0].size(), 1);
//   // EXPECT_EQ(orbits[1].prototype().size(), 1); // 32 (=2*2*2*4) point clusters
//   // EXPECT_EQ(orbits[1].size(), 32);
//   // EXPECT_EQ(orbits[2].prototype().size(), 2); // 6*32 1NN pair clusters
//   // EXPECT_EQ(orbits[2].size(), 6 * 32);
//
//   auto orbit_generators = cluster_specs.make_within_scel_orbit_generators(log());
//
//   jsonParser json;
//   json["orbit_generators"].put_array();
//   for(auto const &cluster_site_indices: orbit_generators) {
//     json["orbit_generators"].push_back(cluster_site_indices);
//   }
//   std::cout << json << std::endl;
//
// }
//
// TEST_F(ClusterSpecsTest, WithinScelTest3) {
//
//   // 3x3x3 conventional FCC (4 atom unit cell / 108 atom supercell) clusters "within_scel"
//   Eigen::Matrix3l T = 3 * _fcc_conventional_transf_mat();
//   Lattice super_lattice = make_superlattice(shared_prim->lattice(), T);
//   SupercellSymInfo supercell_sym_info = make_supercell_sym_info(*shared_prim, super_lattice);
//   auto generating_group = make_generating_group(supercell_sym_info);
//   jsonParser json;
//
//   // a = primitive FCC unit cell lattice vector length
//   double a = shared_prim->lattice()[0].norm();
//
//   std::vector<double> max_length {
//     0, // null cluster
//     0, // point clusters
//     6*a + TOL // 1NN pair clusters only
//   };
//
//   PeriodicMaxLengthClusterSpecs periodic_cluster_specs {
//     shared_prim,
//     shared_prim->factor_group(),
//     alloy_sites_filter,
//     max_length };
//
//   auto periodic_orbits = periodic_cluster_specs.make_periodic_orbits(log());
//   auto scel_periodic_generators = make_scel_periodic_orbit_generators(
//     shared_prim,
//     supercell_sym_info,
//     generating_group,
//     prototype_iterator(periodic_orbits.begin()),
//     prototype_iterator(periodic_orbits.end()));
//   json["scel_periodic_orbit_generators_via_periodic"].put_array();
//   for(auto const &cluster: scel_periodic_generators) {
//     json["orbit_generators"].push_back(make_cluster_site_indices(cluster, supercell_sym_info));
//   }
//
//   auto within_scel_generators_via_periodic = make_within_scel_orbit_generators(
//     supercell_sym_info,
//     generating_group,
//     scel_periodic_generators.begin(),
//     scel_periodic_generators.end());
//
//   json["orbit_generators_via_periodic"].put_array();
//   for(auto const &cluster_site_indices: within_scel_generators_via_periodic) {
//     json["orbit_generators"].push_back(cluster_site_indices);
//   }
//   std::cout << json << std::endl;
//
//   WithinScelMaxLengthClusterSpecs within_scel_cluster_specs {
//     shared_prim,
//     &supercell_sym_info,
//     generating_group,
//     alloy_sites_filter,
//     max_length };
//
//   // auto orbits = cluster_specs.make_within_scel_orbit(log());
//   //
//   // // Construct OrbitPrinterOptions:
//   // OrbitPrinterOptions printer_options;
//   // printer_options.coord_type = INTEGRAL;
//   // // printer_options.print_equivalence_map = true;
//   // // printer_options.print_invariant_group = true;
//   //
//   // // Print prototype clusters and all equivalent clusters
//   // printer_options.orbit_print_mode = ORBIT_PRINT_MODE::PROTO;
//   // print_clust(orbits.begin(), orbits.end(), log(), printer_options);
//   //
//   // EXPECT_EQ(orbits.size(), 3);
//   //
//   // EXPECT_EQ(orbits[0].prototype().size(), 0); // 1 null cluster
//   // EXPECT_EQ(orbits[0].size(), 1);
//   // EXPECT_EQ(orbits[1].prototype().size(), 1); // 108 (=3*3*3*4) point clusters
//   // EXPECT_EQ(orbits[1].size(), 108);
//   // EXPECT_EQ(orbits[2].prototype().size(), 2); // 6*108 1NN pair clusters
//   // EXPECT_EQ(orbits[2].size(), 6 * 108);
//
//   auto orbit_generators = within_scel_cluster_specs.make_within_scel_orbit_generators(log());
//
//   json["orbit_generators"].put_array();
//   for(auto const &cluster_site_indices: orbit_generators) {
//     json["orbit_generators"].push_back(cluster_site_indices);
//   }
//   std::cout << json << std::endl;
//
// }
