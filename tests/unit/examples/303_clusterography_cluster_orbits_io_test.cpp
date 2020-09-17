#include "gtest/gtest.h"

#include "casm/app/AppIO_impl.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"

#include "crystallography/TestStructures.hh" // for test::ZrO_prim

// A test fixture that generates orbits that will be used in the IO examples
class ExampleClusterographyClusterOrbitsIO : public testing::Test {
protected:
  void SetUp() override {

    // Construct a ZrO prim
    auto shared_prim = std::make_shared<CASM::Structure const>(test::ZrO_prim());

    // Use the complete factor group as the generating group
    CASM::SymGroup const &generating_group = shared_prim->factor_group();

    // Cluster "max_length" cutoff values: max_length[b] is cutoff value for orbit branch b
    // - null and point cluster values are ignored
    // - pair clusters up to max length 5.17=c-axis length + a little
    std::vector<double> max_length = {0, 0, 5.17};

    // PeriodicMaxLengthClusterSpecs is constructed with:
    cluster_specs = notstd::make_unique<CASM::PeriodicMaxLengthClusterSpecs>(
                      shared_prim,                        // the prim Structure
                      generating_group,                   // the orbit generating group
                      CASM::alloy_sites_filter,           // filter to include sites with >1 occupant in cluster orbits
                      max_length                          // cluster "max_length" cutoff values
                    );

    orbits = cluster_specs->make_periodic_orbits(CASM::null_log());

    return;
  }

  std::unique_ptr<CASM::PeriodicMaxLengthClusterSpecs> cluster_specs;
  std::vector<CASM::PrimPeriodicIntegralClusterOrbit> orbits;

};

// SymInfoOptions and OrbitPrinterOptions allow controlling some aspects of what is printed and how:
//
// // Controls printing of symmetry operations
// struct SymInfoOptions {
//   COORD_TYPE coord_type = FRAC;   // how symmetry operation coordinates are printed
//   double tol = CASM::TOL;
//   Index prec = 7;
//   bool print_matrix_tau = false;  // print symmetry operations matrix and translation vector
// };
//
// // Controls printing of cluster orbits
// struct OrbitPrinterOptions:
//   int indent_space = 6;
//   char delim = '\n';
//   int prec = 7;
//   COORD_TYPE coord_type = FRAC;
//   ORBIT_PRINT_MODE orbit_print_mode = ORBIT_PRINT_MODE::PROTO;
//   SymInfoOptions sym_info_opt;
//   bool print_coordinates = true;         // print cluster site coordinates
//   bool print_equivalence_map = false;    // print symmetry operations in the equivalence map
//   bool print_invariant_group = false;      // print symmetry operations in the cluster invariant group
// };


TEST_F(ExampleClusterographyClusterOrbitsIO, PrintToStream) {

  // Orbits can be printed to an ostream via CASM::Log.

  // Construct OrbitPrinterOptions:
  CASM::OrbitPrinterOptions printer_options;
  printer_options.coord_type = CASM::INTEGRAL;
  printer_options.print_equivalence_map = true;
  printer_options.print_invariant_group = true;

  // There are two "OrbitPrinter" classes:
  // - ProtoSitesPrinter: print the orbit prototype cluster only
  // - FullSitesPrinter: print all clusters in the orbit

  // The "print_clust" function prints a range of orbits to a stream, using the provided options

  // Print prototype clusters only
  printer_options.orbit_print_mode = CASM::ORBIT_PRINT_MODE::PROTO;
  print_clust(orbits.begin(), orbits.end(), CASM::log(), printer_options);

  // Print prototype clusters and all equivalent clusters
  printer_options.orbit_print_mode = CASM::ORBIT_PRINT_MODE::FULL;
  print_clust(orbits.begin(), orbits.end(), CASM::log(), printer_options);

  EXPECT_EQ(orbits.size(), 6);
}

TEST_F(ExampleClusterographyClusterOrbitsIO, ReadWriteJSON) {

  // Orbits can be written to JSON and read from JSON.

  // Construct OrbitPrinterOptions:
  CASM::OrbitPrinterOptions printer_options;
  printer_options.coord_type = CASM::INTEGRAL;
  printer_options.print_equivalence_map = true;
  printer_options.print_invariant_group = true;

  // The "write_clust" function writes a range of orbits to JSON, using the provided options

  // Print prototype clusters only
  CASM::jsonParser json;
  printer_options.orbit_print_mode = CASM::ORBIT_PRINT_MODE::PROTO;
  write_clust(orbits.begin(), orbits.end(), json, printer_options);

  // Print prototype clusters and all equivalent clusters
  CASM::jsonParser full_orbit_json;
  printer_options.orbit_print_mode = CASM::ORBIT_PRINT_MODE::FULL;
  write_clust(orbits.begin(), orbits.end(), full_orbit_json, printer_options);

  // The JSON object can be written to stream
  CASM::log() << json << std::endl;
  CASM::log() << full_orbit_json << std::endl;

  // To read the orbits back in from JSON we need the prototypes, prim, generating group,
  //   and sym_compare functor
  std::vector<CASM::PrimPeriodicIntegralClusterOrbit> orbits_in;
  read_clust(
    std::back_inserter(orbits_in),
    json,
    *cluster_specs->shared_prim,
    cluster_specs->shared_prim->factor_group(),
    cluster_specs->sym_compare);

  EXPECT_EQ(orbits_in.size(), 6);

  // The full orbit JSON object can be read in similarly (only the prototypes get used)
  std::vector<CASM::PrimPeriodicIntegralClusterOrbit> full_orbit_in;
  read_clust(
    std::back_inserter(full_orbit_in),
    full_orbit_json,
    *cluster_specs->shared_prim,
    cluster_specs->shared_prim->factor_group(),
    cluster_specs->sym_compare);

  EXPECT_EQ(full_orbit_in.size(), 6);
}
