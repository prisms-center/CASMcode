#ifndef CASM_DiffTransConfigEnumPerturbations
#define CASM_DiffTransConfigEnumPerturbations

#include "casm/container/InputEnumerator.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/symmetry/OrbitGeneration.hh"
#include "casm/kinetics/DiffTransEnumEquivalents.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"

namespace CASM {

  namespace Kinetics {

    class Perturbation : public std::set<OccupationTransformation> {

    public:

      Perturbation(std::set<OccupationTransformation> &from_set) {
        for(const OccupationTransformation &item : from_set) {
          this->insert(item);
        }
      };

      Perturbation &apply_sym(const SymOp &op) {
        std::set<OccupationTransformation> tmp;
        for(const OccupationTransformation &occ_trans : *this) {
          tmp.insert(copy_apply(op, occ_trans));
        }
        this->clear();
        for(const OccupationTransformation &item : tmp) {
          this->insert(item);
        }
      };


    };


    /// \brief Enumerate DiffTransConfiguration for a particular DiffusionTransformation,
    ///        set of local clusters, and a particular initial Configuration
    ///
    class DiffTransConfigEnumPerturbations : public InputEnumeratorBase<DiffTransConfiguration> {

      // -- Required members -------------------
      //class Perturbation;
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

      static int run(PrimClex &primclex, const jsonParser &_kwargs, const Completer::EnumOption &enum_opt);

    private:

      /// Implements increment: generate the next DiffTransConfiguration
      void increment() override;

      ///Implements run generates all DiffTransConfigurations given the filter expressions and keyword arguments


      ///... members necessary to do the enumeration ...
      void _init_base_config();
      void _init_unique_difftrans();
      void _init_perturbations();


      void _set_current();
      void _increment_base_config();

      /// Select unique DiffusionTransformations from PrimPeriodicDiffTransOrbit
      //DiffusionTransformation UniqueDiffusionTransformations(const PrimPeriodicDiffTransOrbit &diff_trans_orbit);


      std::set<Perturbation> m_perturbations;
      std::set<Perturbation>::iterator m_perturb_it;

      std::set<DiffusionTransformation>::iterator m_unique_difftrans_it;
      std::set<DiffusionTransformation> m_unique_difftrans;
      DiffTransConfiguration m_base_config;

      Configuration m_background_config;
      PrimPeriodicDiffTransOrbit m_diff_trans_orbit;
      jsonParser m_local_bspecs;
      notstd::cloneable_ptr<DiffTransConfiguration> m_current;
    };

    bool has_local_bubble_overlap(std::vector<LocalIntegralClusterOrbit> &local_orbits, const Supercell &scel);

    std::vector<Supercell> viable_supercells(std::vector<LocalIntegralClusterOrbit> &local_orbits, std::vector<Supercell>);


  }
}
#endif
