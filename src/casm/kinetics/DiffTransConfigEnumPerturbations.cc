#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/clex/Configuration.hh"
#include "casm/kinetics/DiffTransEnumEquivalents.hh"
#include "casm/kinetics/DiffTransConfigEnumPerturbations.hh"
#include "casm/clex/Supercell.hh"

namespace CASM {

  namespace Kinetics {

    DiffTransConfigEnumPerturbations::DiffTransConfigEnumPerturbations(
      const Configuration &background_config,
      const PrimPeriodicDiffTransOrbit &diff_trans_orbit, // or const DiffusionTransformation &diff_trans
      const jsonParser &local_bspecs // or iterators over IntegralClusters
    ) :
      m_background_config(background_config), m_diff_trans_orbit(diff_trans_orbit), m_local_bspecs(local_bspecs), m_dtconfig(background_config, diff_trans_orbit.prototype()) {
    }

    const std::string DiffTransConfigEnumPerturbations::enumerator_name = "DiffTransConfigEnumPerturbations";

    void DiffTransConfigEnumPerturbations::increment() {

      /// Find primitive cell of background_config
      Configuration bg_config_prim = m_background_config.primitive();

      /// Find prototype of m_diff_trans_orbit
      DiffusionTransformation diff_trans_prototype = m_diff_trans_orbit.prototype();
      std::vector<PermuteIterator> permits = bg_config_prim.factor_group();
      std::vector<SymOp> config_fg;
      for(PermuteIterator i : permits) {
        config_fg.push_back(i.sym_op());
      }
      /// Find unique DiffusionTransformations
      //PermuteIterator begin = config_fg.begin();
      //PermuteIterator end = config_fg.end();
      auto begin = bg_config_prim.supercell().permute_begin();
      auto end = bg_config_prim.supercell().permute_end();
      DiffTransEnumEquivalents diff_trans_unique(diff_trans_prototype, begin, end, bg_config_prim);


      /// Use settings from local_bspecs to add perturbations to each from_config
      /// Use DiffTransGroup to check which perturbations are equivalent

    };


  }
}
