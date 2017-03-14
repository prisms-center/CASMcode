#ifndef CASM_DiffTransConfigEnumPerturbations
#define CASM_DiffTransConfigEnumPerturbations

#include "casm/container/InputEnumerator.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/symmetry/OrbitGeneration.hh"
#include "casm/kinetics/DiffTransEnumEquivalents.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"

namespace CASM {

  namespace Kinetics {

    /// \brief Enumerate DiffTransConfiguration for a particular DiffusionTransformation,
    ///        set of local clusters, and a particular initial Configuration
    ///
    class DiffTransConfigEnumPerturbations : public InputEnumeratorBase<DiffTransConfiguration> {

      // -- Required members -------------------

    public:

      /// \brief Construct with an IntegralCluster
      DiffTransConfigEnumPerturbations(
        const Configuration &background_config,
        const PrimPeriodicDiffTransOrbit &diff_trans_orbit, // or const DiffusionTransformation &diff_trans
        const jsonParser &local_bspecs // or iterators over IntegralClusters
      );

      std::string name() const override {
        return enumerator_name;
      }

      static const std::string enumerator_name;

    private:

      /// Implements increment: generate the next DiffTransConfiguration
      void increment() override;

      ///... members necessary to do the enumeration ...

      DiffTransConfiguration m_dtconfig;

      /// Select unique DiffusionTransformations from PrimPeriodicDiffTransOrbit
      //DiffusionTransformation UniqueDiffusionTransformaitions(const PrimPeriodicDiffTransOrbit &diff_trans_orbit);

      Configuration m_background_config;

      PrimPeriodicDiffTransOrbit m_diff_trans_orbit;

      jsonParser m_local_bspecs;

    };

  }
}
#endif
