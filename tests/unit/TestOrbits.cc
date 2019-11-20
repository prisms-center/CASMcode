#include "TestOrbits.hh"
#include <algorithm>
#include <iterator>
#include "casm/app/enum.hh"
#include "casm/kinetics/DiffusionTransformation_impl.hh"
#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"

namespace test {

  TestPrimPeriodicIntegralClusterOrbits::TestPrimPeriodicIntegralClusterOrbits(
    const PrimClex &primclex,
    jsonParser _specs) :
    specs(_specs) {

    try {
      // Make PrimPeriodicIntegralClusterOrbit
      make_prim_periodic_orbits(
        primclex.prim(),
        specs,
        alloy_sites_filter,
        primclex.crystallography_tol(),
        std::back_inserter(orbits),
        primclex.log());
    }
    catch(std::exception &e) {
      default_err_log().error("'make_prim_periodic_orbits' failed in TestPrimPeriodicIntegralClusterOrbits");
      //default_err_log() << "prim: \n" << primclex.prim() << std::endl;
      default_err_log() << "specs: \n" << specs << std::endl;
      default_err_log() << e.what() << std::endl;
      throw e;
    }
  }


  /// Make PrimPeriodicDiffusionTransformationOrbit for all
  ///   PrimPeriodicIntegralClusterOrbit branches
  TestPrimPeriodicDiffusionTransformationOrbits::TestPrimPeriodicDiffusionTransformationOrbits(
    const PrimClex &primclex,
    jsonParser _specs) :
    TestPrimPeriodicIntegralClusterOrbits(primclex, _specs) {

    try {
      // Make PrimPeriodicDiffTransOrbit
      Kinetics::make_prim_periodic_diff_trans_orbits(
        orbits.begin(),
        orbits.end(),
        primclex.crystallography_tol(),
        std::back_inserter(diff_trans_orbits),
        &primclex);
    }
    catch(std::exception &e) {
      default_err_log().error("'Kinetics::make_prim_periodic_diff_trans_orbits' failed in TestPrimPeriodicDiffusionTransformationOrbits");
      //default_err_log() << "prim: \n" << primclex.prim() << std::endl;
      default_err_log() << "specs: \n" << specs << std::endl;
      default_err_log() << "orbits.size(): " << orbits.size() << std::endl;
      default_err_log() << e.what() << std::endl;
      throw e;
    }
  }

  /// Make PrimPeriodicDiffusionTransformationOrbit for a range of
  ///   PrimPeriodicIntegralClusterOrbit branches
  TestPrimPeriodicDiffusionTransformationOrbits::TestPrimPeriodicDiffusionTransformationOrbits(
    const PrimClex &primclex,
    jsonParser _specs,
    Index branch_begin,
    Index branch_end) :
    TestPrimPeriodicIntegralClusterOrbits(primclex, _specs) {

    try {

      auto begin = orbits.begin();
      for(auto it = begin; it != orbits.end(); ++it) {
        if(it->prototype().size() <= branch_begin) {
          begin = it;
          break;
        }
      }

      auto end = orbits.end();
      if(begin != end) {
        auto it = begin;
        ++it;
        for(; it != orbits.end(); ++it) {
          if(it->prototype().size() >= branch_end) {
            end = it;
            break;
          }
        }
      }

      // Make PrimPeriodicDiffTransOrbit
      Kinetics::make_prim_periodic_diff_trans_orbits(
        begin,
        end,
        primclex.crystallography_tol(),
        std::back_inserter(diff_trans_orbits),
        &primclex);
    }
    catch(std::exception &e) {
      default_err_log().error("'Kinetics::make_prim_periodic_diff_trans_orbits' failed in TestPrimPeriodicDiffusionTransformationOrbits");
      //default_err_log() << "prim: \n" << primclex.prim() << std::endl;
      default_err_log() << "specs: \n" << specs << std::endl;
      default_err_log() << "orbits.size(): " << orbits.size() << std::endl;
      default_err_log() << "branch_begin: " << branch_begin << std::endl;
      default_err_log() << "branch_end: " << branch_end << std::endl;
      default_err_log() << e.what() << std::endl;
      throw e;
    }
  }

  /// Make PrimPeriodicDiffusionTransformationOrbit using DiffusionTransformation::run
  TestPrimPeriodicDiffusionTransformationOrbits::TestPrimPeriodicDiffusionTransformationOrbits(
    const PrimClex &primclex,
    jsonParser diff_trans_specs,
    const Completer::EnumOption &enum_opt) :
    TestPrimPeriodicIntegralClusterOrbits(primclex, diff_trans_specs["cspecs"]) {

    try {
      // Generate DiffTrans
      std::set<PrimPeriodicDiffTransOrbit> tmp;
      Kinetics::DiffusionTransformationEnum::run(primclex, diff_trans_specs, enum_opt, tmp);
      std::copy(tmp.begin(), tmp.end(), std::inserter(diff_trans_orbits, diff_trans_orbits.end()));
    }
    catch(std::exception &e) {
      default_err_log().error("'Kinetics::DiffusionTransformationEnum::run' failed in TestPrimPeriodicDiffusionTransformationOrbits");
      //default_err_log() << "prim: \n" << primclex.prim() << std::endl;
      default_err_log() << "diff_trans_specs: \n" << diff_trans_specs << std::endl;
      default_err_log() << e.what() << std::endl;
      throw e;
    }
  }


  TestLocalOrbits::TestLocalOrbits(
    const PrimClex &primclex,
    const Kinetics::DiffusionTransformation &diff_trans,
    jsonParser _specs) :
    specs(_specs) {

    try {
      SymGroup generating_grp {
        diff_trans.invariant_subgroup(
          primclex.prim().factor_group(),
          PrimPeriodicDiffTransSymCompare(primclex.shared_prim(), primclex.crystallography_tol()))};
      LocalSymCompare<IntegralCluster> sym_compare(primclex.shared_prim(), primclex.crystallography_tol());

      make_local_orbits(
        diff_trans,
        generating_grp,
        sym_compare,
        specs,
        alloy_sites_filter,
        primclex.crystallography_tol(),
        std::back_inserter(orbits),
        primclex.log());
    }
    catch(std::exception &e) {
      default_err_log().error("'make_local_orbits' failed in TestLocalOrbits");
      //default_err_log() << "prim: \n" << primclex.prim() << std::endl;
      default_err_log() << "diff_trans: \n" << diff_trans << std::endl;
      default_err_log() << "specs: \n" << specs << std::endl;
      default_err_log() << e.what() << std::endl;
      throw e;
    }
  }
}
