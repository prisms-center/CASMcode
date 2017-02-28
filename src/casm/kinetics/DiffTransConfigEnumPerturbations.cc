#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {

  namespace Kinetics {

  DiffTransConfigEnumPerturbations(
        const Configuration &background_config,
        const PrimPeriodicDiffTransOrbit &diff_trans_orbit, // or const DiffusionTransformation &diff_trans
        const jsonParser &local_bspecs // or iterators over IntegralClusters
      ) :
        m_background_config(background_config), m_diff_trans_orbit(diff_trans_orbit), m_local_bspecs(local_bspecs){
        }





  /// Find primitive cell of background_config
  Configuration prim_bg_config = background.primitive();

  /// Determine factor group of prim_bg_config
  std::vector<PermuteIterator> f_group = prim_bg_config.factor_group();

  /// Find all possible translations of diff_trans in the supercell of the primitive config
  DiffusionTransformation translations = prim_bg_config.supercell().prim_grid().unitcell();

  /// Iterate over DiffusionTransformation in diff_trans_orbit
  for (auto it = diff_trans_orbit.begin(); it != diff_trans_orbit.end(); ++it)
  {
    DiffusionTransformation diff_trans = diff_trans_orbit[it];

    /// Apply translations to each diff_trans in diff_trans_orbit
    for (auto it2 = translations.begin(); it2 != translations.end(); ++it2)
    {

      /// Apply f_group to each translation

        /// If all symmetry operations from f_group result in unique DiffusionTransformations, then keep it.

        /// If any of the symmetry operations map the translation onto an existing one,
        /// translate the DiffusionTransformation back into the supercell and then compare to determine if it is unique
        /// ^ Have to determine how to do this in a smart way. Anton suggested using clusters within a certain radius

        /// If diffusion transformations are still non-unique, compare them and keep the greater one.
        /// Use Configuration factor group & ScelPeriodicSymCompare symmetry
        /// ScelPeriodicSymCompare< IntegralCluster >::ScelPeriodicSymCompare(const PrimGrid &prim_grid, double tol)

    }
  }
  /// Make from_config for each unique DiffusionTransformation + background_config


  /// For each unique DiffusionTransformation, find the DiffTransGroup, which is the invariant subgroup of Configuration::factor_group()
  /// that leaves the ‘from_config’ unchanged (both the DiffusionTransformation and the rest of the occupants unchanged)


  /// Use settings from local_bspecs to add perturbations to each from_config
  /// Use DiffTransGroup to check which perturbations are equivalent

  }
}
