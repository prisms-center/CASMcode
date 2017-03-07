#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/clex/Configuration.hh"
#include "casm/kinetics/DiffTransEnumEquivalents.hh"

namespace CASM {

  namespace Kinetics {

    DiffTransConfigEnumPerturbations(
      const Configuration &background_config,
      const PrimPeriodicDiffTransOrbit &diff_trans_orbit, // or const DiffusionTransformation &diff_trans
      const jsonParser &local_bspecs // or iterators over IntegralClusters
    ) :
      m_background_config(background_config), m_diff_trans_orbit(diff_trans_orbit), m_local_bspecs(local_bspecs) {
    }





    /// Find primitive cell of background_config
    Configuration bg_config_prim = m_background_config.primitive();

    /// Find prototype of m_diff_trans_orbit
    DiffusionTransformation diff_trans_prototype = diff_trans_orbit.prototype();

    /// Find unique DiffusionTransformations

    DiffTransEnumEquivalents diff_trans_unique(diff_trans_prototype, begin, end, bg_config_prim);


    /// Use settings from local_bspecs to add perturbations to each from_config
    /// Use DiffTransGroup to check which perturbations are equivalent

  }
}
