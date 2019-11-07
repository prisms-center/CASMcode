#include "gtest/gtest.h"
#include "autotools.hh"

/// What is being tested:
#include "casm/kinetics/DiffTransConfigEnumOccPerturbations.hh"

/// What is being used to test it:
#include "casm/clex/PrimClex.hh"
#include "Common.hh"
#include "TestConfiguration.hh"
#include "TestOrbits.hh"
#include "casm/app/AppIO_impl.hh"
#include "casm/app/enum.hh"
#include "casm/app/QueryHandler.hh"
#include "casm/casm_io/jsonFile.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/ConfigEnumAllOccupations_impl.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/container/multivector.hh"
#include "casm/kinetics/DiffusionTransformation_impl.hh"
#include "casm/kinetics/DiffTransConfiguration_impl.hh"
#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"

using namespace CASM;
using namespace test;

typedef Orbit <
//Kinetics::DiffusionTransformation,
Kinetics::PrimPeriodicDiffTransSymCompare > PrimPeriodicDiffTransOrbit;

namespace {

  struct TestOrbits0 : test::TestPrimPeriodicDiffusionTransformationOrbits {
    TestOrbits0(const PrimClex &primclex) :
      test::TestPrimPeriodicDiffusionTransformationOrbits(
        primclex,
        jsonFile(autotools::abs_srcdir() + "/tests/unit/kinetics/ZrO_bspecs_0.json"),
        2, 4) { // orbit_branches [2,4)
      EXPECT_EQ(orbits.size(), 74);
      EXPECT_EQ(diff_trans_orbits.size(), 4);
    }
  };

  struct TestUnitConfig0 : test::TestConfiguration {

    TestUnitConfig0(const PrimClex &primclex) :
      TestConfiguration(
        primclex,
        Eigen::Vector3i(1, 1, 1).asDiagonal(), {
      0, 0, 1, 0
    }) {}

  };

  struct TestBackgroundConfig0 : test::TestConfiguration {

    TestBackgroundConfig0(const PrimClex &primclex) :
      TestConfiguration(primclex, background_config(primclex)) {}

    static Configuration background_config(const PrimClex &primclex) {
      Supercell background_scel {&primclex, Eigen::Matrix3i(Eigen::Vector3i(3, 3, 3).asDiagonal())};
      TestUnitConfig0 tu(primclex);
      return tu.config.fill_supercell(background_scel, primclex.prim().factor_group()).in_canonical_supercell();
    }

  };

  struct TestEnumComponents {

    typedef Kinetics::DiffTransConfigEnumOccPerturbations EnumType;
    typedef EnumType::CurrentPerturbation CurrentPerturbation;

    Log &log;
    std::unique_ptr<EnumType> enumerator;

    // for each base, collect multiplicity of local orbits:
    //   local_orbit_size[base_index][local_orbit_index] =
    //     enumerator->local_orbit()[local_orbit_index].size()
    multivector<int>::X<2> local_orbit_mult;
    multivector<int>::X<2> local_orbit_branch;

    // for each (base, orbit, occ), collect proposed perturbation:
    //   all_perturb[base_index][local_orbit_index][occ_counter_index] =
    //      enumerator->current_perturb()
    multivector<CurrentPerturbation>::X<3> all_perturb; //base, orbit, occ

    // for each base, collect unique canonical perturbations:
    //    unique_canon_perturb[base] = enumerator->current_perturb().perturb
    multivector<OccPerturbation>::X<2> unique_canon_perturb;

    TestEnumComponents(
      const Configuration &bg_config,
      const PrimPeriodicDiffTransOrbit &dt_orbit,
      jsonParser local_specs,
      Log &_log = null_log()):
      log(_log) {

      log.begin("TestEnumComponents");
      log << "bg_config: \n" << bg_config << std::endl;
      log << "diff_trans: \n" << dt_orbit.prototype() << std::endl;
      log << "local_specs: \n" << local_specs << std::endl;

      /// Construct enumerator
      log.construct("Kinetics::DiffTransConfigEnumOccPerturbations");
      enumerator.reset(new Kinetics::DiffTransConfigEnumOccPerturbations(
                         bg_config,
                         dt_orbit,
                         local_specs));
      log << "done" << std::endl;
      EXPECT_EQ(true, true);
      EXPECT_EQ(enumerator->valid(), true);
      EXPECT_EQ(enumerator->step(), 0);
      EXPECT_EQ((enumerator->begin() != enumerator->end()), true);
      log << "basic checks ok" << std::endl;

      run();
    }

    void print_info() {
      log << "base: " << enumerator->base_index()
          << "/" << enumerator->base().size()
          << "  local_orbit: " << enumerator->local_orbit_index()
          << "/" << enumerator->local_orbit().size()
          << "  occupation: " << enumerator->occ_counter_index()
          << std::endl;
      const auto &curr = enumerator->current_perturb();
      log << "  current_perturb: "
          << curr.is_not_subcluster << " " << curr.is_canonical
          << " :: valid perturbation?: " << (curr.is_not_subcluster && curr.is_canonical) << std::endl
          << curr.perturb;
    }

    bool step() {

      const auto &base = enumerator->base();
      Index base_i = enumerator->base_index();
      Index orbit_i = enumerator->local_orbit_index();
      Index occ_i = enumerator->occ_counter_index();

      if(enumerator->valid()) {
        // first time at base, print info
        if(base_i == all_perturb.size()) {
          log << "base[" << base_i << "].config: \n"
              << base[base_i].config << std::endl;
          log << "base[" << base_i << "].diff_trans: \n"
              << base[base_i].diff_trans << std::endl;
          log << "base[" << base_i << "].diff_tran_g.size(): "
              << base[base_i].diff_trans_g.size() << std::endl;
          log << "base[" << base_i << "].generating_sym_g.size(): "
              << base[base_i].generating_sym_g.size() << std::endl;

          // print symgroup
          const auto &g = base[base_i].generating_sym_g;
          auto frac = [&](const Eigen::VectorXd & vec) {
            return base[base_i].diff_trans.prim().lattice().lat_column_mat().inverse() * vec;
          };
          for(Index op_i = 0; op_i < g.size(); ++op_i) {
            log << "  op " << op_i << "/" << g.size() << ":\n"
                << "master group index: " << g[op_i].index() << std::endl
                << "matrix: \n" << g[op_i].matrix() << std::endl
                << "tau: " << frac(g[op_i].tau()).transpose() << std::endl
                << "integral_tau: " << frac(g[op_i].integral_tau()).transpose()
                << std::endl << std::endl;
          }
        }

        // first time at local orbit, print clusters
        if(occ_i == 0) {
          log << "local_orbit " << orbit_i << "/" << enumerator->local_orbit().size() << ":\n";
          Index clust_i = 0;
          for(const auto &clust : enumerator->local_orbit()[orbit_i]) {
            log << "cluster " << clust_i++ << "/" << enumerator->local_orbit()[orbit_i].size() << ":\n"
                << clust;
          }
          log << std::endl;

          log << "check " << orbit_i << "/" << enumerator->local_orbit().size() << ":\n";
          const auto &g = base[base_i].generating_sym_g;
          auto proto = enumerator->local_orbit()[orbit_i].prototype();
          ScelPeriodicSymCompare<IntegralCluster> sym_compare(base[base_i].config.supercell());

          for(Index op_i = 0; op_i < g.size(); ++op_i) {
            log << "op " << op_i << "/" << g.size() << ":\n"
                << sym_compare.prepare(copy_apply(g[op_i], proto));
          }
          log << std::endl;
        }
      }

      print_info();

      log << "do check increment" << std::endl;
      bool check = enumerator->check_increment();
      log << "  done: " << check << std::endl;

      if(enumerator->valid()) {

        // expand all_perturb as necessary and add current perturb
        if(base_i == all_perturb.size()) {
          all_perturb.push_back(multivector<CurrentPerturbation>::X<2>());
        }
        if(orbit_i == all_perturb.back().size()) {
          all_perturb.back().push_back(multivector<CurrentPerturbation>::X<1>());
        }
        all_perturb[base_i][orbit_i].push_back(enumerator->current_perturb());

        // collect multiplicity, branch of local orbits, for each base config
        if(base_i == local_orbit_mult.size()) {
          local_orbit_mult.push_back(multivector<int>::X<1>());
          local_orbit_branch.push_back(multivector<int>::X<1>());
          for(const auto &local_orbit : enumerator->local_orbit()) {
            local_orbit_mult[base_i].push_back(local_orbit.size());
            local_orbit_branch[base_i].push_back(local_orbit.prototype().size());
          }
        }

        // collect unique, canonical perturbations, for each base config
        if(check) {
          if(base_i == unique_canon_perturb.size()) {
            unique_canon_perturb.push_back(multivector<OccPerturbation>::X<1>());
          }
          unique_canon_perturb.back().push_back(enumerator->current_perturb().perturb);
        }
      }

      return check;
    }

    void run() {
      log.custom("Initialization");
      log << std::endl << std::endl;

      step();

      do {
        log.begin("Increment");
        log << std::endl;

        bool check;
        do {

          log.begin("Partial Increment");
          log << "do partial increment" << std::endl;
          enumerator->partial_increment(true);
          log << "  done" << std::endl;

          check = step();

          log << std::endl;

        }
        while(!check);

        log.custom("Finish Increment");
        log << "enumerator valid: " << enumerator->valid() << std::endl << std::endl;

      }
      while(enumerator->valid());

      log.end("TestEnumComponents");
      EXPECT_EQ(enumerator->valid(), false);

    }
  };

}


TEST(DiffTransConfigEnumOccPerturbationsTest, NeighborhoodOverlapTest) {

  /// Make test project
  EXPECT_EQ(true, true);
  test::ZrOProj proj;
  proj.check_init();
  proj.check_composition();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);
  //const Structure &prim = primclex.prim();
  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();
  EXPECT_EQ(true, true);

  // Make orbits
  TestOrbits0 to(primclex);

  // Make local orbits
  TestLocalOrbits local_to(
    primclex,
    to.diff_trans_orbits[0].prototype(),
    jsonFile(autotools::abs_srcdir() + "/tests/unit/kinetics/ZrO_local_bspecs_1.json"));

  ///Make various supercells
  Supercell scel1 {&primclex, Lattice(2 * a, 2 * b, 3 * c)};
  Supercell scel2 {&primclex, Lattice(4 * a, 4 * b, 3 * c)};
  Supercell scel3 {&primclex, Lattice(8 * a, 2 * b, 3 * c)};
  Supercell scel4 {&primclex, Lattice(4 * a, 3 * b, 4 * c)};
  std::vector<Supercell> scel_list;
  scel_list.push_back(scel1);
  scel_list.push_back(scel2);
  scel_list.push_back(scel3);
  scel_list.push_back(scel4);

  EXPECT_EQ(has_local_neighborhood_overlap(local_to.orbits, scel1), 1);
  EXPECT_EQ(has_local_neighborhood_overlap(local_to.orbits, scel2), 0);
  EXPECT_EQ(has_local_neighborhood_overlap(local_to.orbits, scel3), 1);
  EXPECT_EQ(has_local_neighborhood_overlap(local_to.orbits, scel4), 1);
  std::vector<Supercell> result = viable_supercells(local_to.orbits, scel_list);
  EXPECT_EQ(*(result.begin()) == scel2, 1);

}

TEST(DiffTransConfigEnumOccPerturbationsTest, ZrOTest_Components) {

  /// Make test project
  test::ZrOProj proj;
  proj.check_init();
  Log &log = null_log();
  PrimClex primclex(proj.dir, null_log());
  EXPECT_EQ(true, true);

  // Make orbits & background config
  log.construct("TestOrbits0");
  TestOrbits0 to(primclex);
  log << "  DONE" << std::endl;

  /// background config 0
  {
    log.construct("TestBackgroundConfig0");
    log << "ZrO 3x3x3 supercell, alternating filled and empty O layers" << std::endl;
    TestBackgroundConfig0 tc(primclex);
    EXPECT_EQ(true, true);
    log << "  DONE" << std::endl;

    /// Step-by-step checks
    {
      /// Common local cluster specs
      jsonFile local_specs(autotools::abs_srcdir() + "/tests/unit/kinetics/ZrO_local_bspecs_0.json");

      {
        log.check("DiffTrans orbit 0");
        log << "hop along c-axis" << std::endl;
        TestEnumComponents te(tc.config, to.diff_trans_orbits[0], local_specs, log);
        log << "- expected base configurations:" << 1 << std::endl;
        log << "  - base 0: c-axis hop between neighboring filled and empty layer" << std::endl;
        EXPECT_EQ(te.enumerator->base().size(), 1);
        log << "- expected local orbit branch, multiplicity:" << std::endl;
        log << "  - base config 0:" << std::endl;
        log << "    - 1: null perturbation" << std::endl;
        log << "    - 6: nearest point cluster in O layer" << std::endl;
        log << "    - 6: nearest point cluster in Va layer" << std::endl;
        EXPECT_EQ(te.local_orbit_branch[0], std::vector<int>({0, 1, 1}));
        EXPECT_EQ(te.local_orbit_mult[0], std::vector<int>({1, 6, 6}));
        log << "- expected unique, canonical DiffTransConfigurations:" << std::endl;
        log << "  - base config 0:" << std::endl;
        log << "    - 0: null perturbation" << std::endl;
        log << "    - 1: nearest Va in O layer" << std::endl;
        log << "    - 2: nearest O in Va layer" << std::endl;
        EXPECT_EQ(te.unique_canon_perturb[0].size(), 3);
        log << "done" << std::endl;

      }
    } // end Step-by-step checks

  }
}

TEST(DiffTransConfigEnumOccPerturbationsTest, ZrOTest_run) {

  /// Make test project
  EXPECT_EQ(true, true);
  test::ZrOProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  //const Structure &prim(primclex.prim());
  primclex.settings().set_crystallography_tol(1e-5);

  // -- Generate Supercell & Configuration --
  Completer::EnumOption enum_opt;
  enum_opt.desc();

  jsonParser kwargs;
  to_json(ScelEnumProps(1, 5), kwargs["supercells"]);
  kwargs["supercells"]["existing_only"] = false;
  //EXPECT_EQ(true, true);

  ConfigEnumAllOccupations::run(primclex, kwargs, enum_opt, nullptr);
  //EXPECT_EQ(true, true);

  // Test Kinetics::DiffTransConfigEnumOccPerturbations::run
  {
    Completer::EnumOption enum_opt;
    enum_opt.desc();

    // Generate DiffTrans
    jsonFile diff_trans_json {autotools::abs_srcdir() + "/tests/unit/kinetics/ZrO_diff_trans_0.json"};
    int success = Kinetics::DiffusionTransformationEnum::run(primclex, diff_trans_json, enum_opt, nullptr);
    const auto &diff_trans_db = primclex.generic_db<Kinetics::PrimPeriodicDiffTransOrbit>();
    EXPECT_EQ(diff_trans_db.size(), 3);
    EXPECT_EQ(success, 0);

    //    //print DiffTrans prototypes
    //    {
    //      PrototypePrinter<Kinetics::DiffusionTransformation> printer;
    //      print_clust(diff_trans_db.begin(), diff_trans_db.end(), std::cout, printer);
    //    }

    // Generate perturbations
    jsonFile diff_perturb_json {autotools::abs_srcdir() + "/tests/unit/kinetics/ZrO_diff_perturb_0.json"};
    Kinetics::DiffTransConfigEnumOccPerturbations::run(primclex, diff_perturb_json, enum_opt, nullptr);
    EXPECT_EQ(true, true);

    // Test (quantity 1856 not checked for accuracy)
    EXPECT_EQ(primclex.generic_db<Kinetics::DiffTransConfiguration>().size(), 1856);
  }
}

