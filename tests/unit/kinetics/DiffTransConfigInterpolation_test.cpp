#include "gtest/gtest.h"

/// What is being tested:
#include "casm/kinetics/DiffTransConfigInterpolation.hh"

/// What is being used to test it:
#include "casm/clex/PrimClex.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/AppIO_impl.hh"
#include "Common.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/kinetics/DoFTransformation.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/DiffusionTransformationEnum.hh"
#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"

using namespace CASM;
using namespace test;


TEST(DiffTransConfigInterpolationTest, Test0) {
  /*
  test::ZrOProj proj;
  proj.check_init();
  proj.check_composition();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);

  fs::path bspecs_path = "tests/unit/kinetics/bspecs_0.json";
  jsonParser bspecs {bspecs_path};

  std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
  make_prim_periodic_orbits(
    primclex.prim(),
    bspecs,
    alloy_sites_filter,
    primclex.crystallography_tol(),
    std::back_inserter(orbits),
    primclex.log());

  //print_clust(orbits.begin(), orbits.end(), primclex.log(), ProtoSitesPrinter());
  std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
  Kinetics::make_prim_periodic_diff_trans_orbits(
    orbits.begin() + 4,
    orbits.begin() + 7,
    primclex.crystallography_tol(),
    std::back_inserter(diff_trans_orbits),
    &primclex);
  Kinetics::DiffusionTransformation trans = diff_trans_orbits[0].prototype();
  Kinetics::DiffusionTransformation trans2 = diff_trans_orbits[2].prototype();

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  Supercell scel {&primclex, Lattice(2 * a, 2 * b, 2 * c)};
  Supercell scel_prim {&primclex, Lattice(a, b, c)};

  // include displacements & deformation
  Eigen::Vector3d zero(0., 0., 0.);
  Eigen::Vector3d dx(0.001, 0., 0.);
  Eigen::Vector3d dy(0., 0.001, 0.);
  Eigen::Vector3d dz(0., 0., 0.001);

  Eigen::Matrix3d dF1;
  dF1 << 0.01, 0.01, 0.0,
    0.0, -0.01, 0.0,
    0.0, 0.0, 0.05;
  Eigen::Matrix3d dF2;
  dF2 << 0.0, 0.02, 0.01,
    -0.01, -0.01, 0.0,
    0.0, -0.01, 0.02;
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();

  Configuration config(scel);
  config.init_occupation();
  config.init_displacement();
  config.init_deformation();
  config.init_specie_id();
  //hardcoded occupation for trans to occur is there a way to do this generally?
  config.set_occupation({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1});
  config.set_disp(5, dx*2);
  config.set_disp(6, dy*-2);
  config.set_disp(27, dz*3);
  config.set_disp(30, dy*1 + dx*4);
  config.set_disp(31, dz*4 + dx*3);

  Configuration config_prim(scel_prim);
  config_prim.init_occupation();
  config_prim.init_specie_id();
  config_prim.set_occupation({0,0,1,1});
  //test Constructor/field accessors
  Kinetics::DiffTransConfiguration dtc(config, trans);

  Configuration to_config = dtc.to_config();

  Kinetics::DiffTransConfigInterpolation dtc_interpol(dtc,4);
  // int temp = dtc_inter.run("abc",4,dtc);

  config.write_pos(std::cout);
  config_prim.write_pos(std::cout);
  // Configuration::displacement_matrix_t config_disp = config.displacement();
  // std::cout << config_disp << "\n";
  // Configuration::displacement_matrix_t to_config_disp = to_config.displacement();
  // std::cout << to_config_disp << "\n";
  */
}

