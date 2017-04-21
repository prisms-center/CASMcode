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
      static const std::string interface_help;

      static int run(PrimClex &primclex, const jsonParser &kwargs, const Completer::EnumOption &enum_opt);

    private:

      /// Implements increment: generate the next DiffTransConfiguration
      void increment() override;

      ///Implements run generates all DiffTransConfigurations given the filter expressions and keyword arguments


      ///... members necessary to do the enumeration ...
      void _init_base_dtc();
      void _init_unique_dts();
      void _init_perturbations();


      void _set_current();
      void _increment_base_dtc();

      /// Select unique DiffusionTransformations from PrimPeriodicDiffTransOrbit
      //DiffusionTransformation UniqueDiffusionTransformaitions(const PrimPeriodicDiffTransOrbit &diff_trans_orbit);
      /*
      more stuff about perturbation index and such
      std::set<PerturbationType> m_perturbations;
      std::set<PerturbationType>::iterator m_perturb_it;
      */

      std::set<DiffusionTransformation>::iterator m_unique_dts_it;
      std::set<DiffusionTransformation> m_unique_dts;
      DiffTransConfiguration m_base_dtc;

      Configuration m_background_config;
      PrimPeriodicDiffTransOrbit m_diff_trans_orbit;
      jsonParser m_local_bspecs;
      notstd::cloneable_ptr<DiffTransConfiguration> m_current;
    };

  }
}
#endif
