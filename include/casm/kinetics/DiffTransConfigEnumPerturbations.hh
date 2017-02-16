#ifndef CASM_DiffTransConfigEnumPerturbations
#define CASM_DiffTransConfigEnumPerturbations

#include "casm/container/InputEnumerator.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/symmetry/OrbitGeneration.hh"

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
        const DiffusionTransformation &diff_trans,  // or const PrimPeriodicDiffTransOrbit &diff_trans_orbit,
        const jsonParser &local_bspecs); // or iterators over IntegralClusters);

      std::string name() const override {
        return enumerator_name;
      }

      static const std::string enumerator_name;

      private:


      /// Implements increment: generate the next DiffTransConfiguration
      void increment() override;

      ... members necessary to do the enumeration ...

        DiffTransConfiguration m_dtconfig;
      };

  }

#endif